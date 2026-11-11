"""Downloadable real docking samples. No synthetic padding of pose counts."""
import csv, io, json, hashlib, math, re, zipfile
from contextvars import ContextVar
from pathlib import Path
from uuid import uuid4
from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, Response
from .config import DATA_DIR

sampling_seed = ContextVar('sampling_seed', default=42)
router = APIRouter()
EXPORTS = DATA_DIR / 'docking_exports'

def pdbqt_records(text):
    blocks = re.findall(r'^MODEL\b.*?^ENDMDL[^\n]*', text, re.M | re.S)
    return blocks or ([text] if text.strip() else [])


def _pose_as_pdb(item):
    """Normalize one engine pose to PDB coordinates for the viewer and RMSD service."""
    if item.get('common_pose_pdb'):
        return str(item['common_pose_pdb'])
    if item.get('pdb'):
        return str(item['pdb'])
    text = str(item.get('text') or '').strip()
    if not text:
        return ''
    fmt = str(item.get('format') or '').lower()
    if fmt == 'pdb':
        return text
    if fmt == 'pdbqt':
        try:
            from .engines.poses import pose_to_pdb
            return pose_to_pdb(text)
        except Exception:
            return ''
    if fmt == 'sdf':
        try:
            from rdkit import Chem
            mol = Chem.MolFromMolBlock(text.split('$$$$', 1)[0], removeHs=True, sanitize=False)
            return Chem.MolToPDBBlock(mol) if mol is not None else ''
        except Exception:
            return ''
    return ''


def dock_with_download(*args, **kwargs):
    from .docking import dock_smiles
    requested_poses = max(1, min(20, int(kwargs.get('n_poses', 10) or 10)))
    records, results, warnings = [], [], []
    for attempt in range(requested_poses):
        token = sampling_seed.set(42 + attempt * 7919)
        try:
            result = dock_smiles(*args, **kwargs)
        except Exception as exc:
            if not results:
                raise
            warnings.append(f'Additional sampling stopped: {type(exc).__name__}: {exc}')
            break
        finally:
            sampling_seed.reset(token)
        results.append(result)
        metadata = result.get('engine_metadata') or {}
        native = metadata.pop('export_poses', None)
        if not native:
            # Vina returns every mode under metadata.poses. Export those modes
            # directly instead of retaining only result.docked_pdb (the best).
            native = [
                {
                    'text': pose.get('pdb') or pose.get('text') or '',
                    'format': pose.get('format') or ('pdb' if pose.get('pdb') else 'pdbqt'),
                    'energy': pose.get('energy'),
                    'score_type': pose.get('score_type') or metadata.get('energy_source') or 'Binding energy (kcal/mol)',
                }
                for pose in (metadata.get('poses') or [])
            ]
        if not native:
            native = [{'text': result.get('docked_pdb') or result.get('docked_pdbqt') or '',
                       'format': 'pdb' if result.get('docked_pdb') else 'pdbqt',
                       'energy': result.get('energy'),
                       'score_type': metadata.get('energy_source') or 'Binding energy (kcal/mol)'}]
        for item in native:
            text = str(item.get('text') or item.get('pdb') or '')
            if not text.strip():
                continue
            record = {**item, 'text': text, 'attempt': attempt + 1,
                      'execution_mode': metadata.get('mode') or ''}
            details = record.get('native_details') or {}
            details = dict(details) if isinstance(details, dict) else {}
            existing_score_names = {str(name).lower().replace('_', '') for name in details}
            for key, value in record.items():
                if value is None or key in {'native_score', 'native_score_type', 'native_details', 'score_type', 'energy', 'common_be'}:
                    continue
                name = str(key).lower()
                if name.replace('_', '') not in existing_score_names and any(token in name for token in ('score', 'affinity', 'energy')):
                    try:
                        numeric = float(value)
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(numeric):
                        details.setdefault(key, numeric)
            record['native_details'] = details
            # Count only poses that can actually be normalized for the viewer
            # and RMSD service; malformed output must not satisfy the requested
            # pose count and suppress another sampling attempt.
            if _pose_as_pdb(record).strip():
                records.append(record)
        if len(records) >= requested_poses:
            break
    if not records:
        raise ValueError('Engine returned no downloadable coordinates.')
    _, best = min(enumerate(results), key=lambda pair: pair[1]['energy'])
    metadata = best.setdefault('engine_metadata', {})
    ensemble = []
    for index, item in enumerate(records, 1):
        pdb = _pose_as_pdb(item)
        if not pdb.strip():
            continue
        energy = item.get('common_be')
        score_type = str(
            (item.get('common_be_source') or item.get('score_type') or 'Docking score')
            if energy is not None else (item.get('score_type') or 'Docking score')
        )
        if energy is None and 'confidence' not in score_type.lower() and 'karma score' not in score_type.lower():
            energy = item.get('energy')
        try:
            energy = float(energy) if energy is not None else None
        except (TypeError, ValueError):
            energy = None
        is_confidence = 'confidence' in score_type.lower()
        ensemble.append({
            'pdb': pdb, 'energy': energy, 'score_type': score_type,
            'pose_index': index, 'confidence': item.get('confidence', item.get('energy')) if is_confidence else item.get('confidence'),
            'native_score': item.get('native_score', item.get('energy')),
            'native_score_type': item.get('native_score_type') or item.get('score_type') or score_type,
            'native_details': item.get('native_details') or {},
            'execution_mode': item.get('execution_mode') or metadata.get('mode') or '',
            'attempt': item.get('attempt'), 'format': item.get('format') or 'pdb',
        })

    best_pose = str(best.get('docked_pdb') or '')
    if best_pose.strip() and not any(p['pdb'].strip() == best_pose.strip() for p in ensemble):
        ensemble.append({
            'pdb': best_pose, 'energy': float(best['energy']),
            'score_type': metadata.get('energy_source') or 'Binding energy (kcal/mol)',
            'pose_index': len(ensemble) + 1, 'confidence': None,
            'native_score': float(best['energy']),
            'native_score_type': metadata.get('energy_source') or 'Binding energy (kcal/mol)',
            'native_details': {}, 'attempt': None, 'format': 'pdb',
        })
    ensemble.sort(key=lambda p: (p['energy'] is None, p['energy'] if p['energy'] is not None else float('inf')))
    for index, pose in enumerate(ensemble, 1):
        pose['pose_index'] = index
    # Use the lowest returned binding-energy pose as the reference. Engines
    # that expose confidence only keep the RMSD blank instead of mislabelling it.
    energy_poses = [pose for pose in ensemble if pose.get('energy') is not None and 'confidence' not in str(pose.get('score_type', '')).lower()]
    if energy_poses:
        try:
            from .poses import pose_mol, rmsd_common
            reference = min(energy_poses, key=lambda pose: float(pose['energy']))
            reference_mol = pose_mol(reference['pdb'])
            reference['rmsd_to_best_pose'] = 0.0
            reference['atoms_compared'] = None
            for pose in ensemble:
                if pose is reference:
                    continue
                comparison = rmsd_common(reference_mol, pose_mol(pose['pdb']), min_atoms=1)
                pose['rmsd_to_best_pose'] = round(float(comparison[0]), 3) if comparison else None
                pose['atoms_compared'] = int(comparison[1]) if comparison else None
        except Exception:
            pass
    for pose in ensemble:
        pose['file'] = f"poses/pose_{int(pose['pose_index']):03d}_normalized.pdb"
    metadata['pose_ensemble'] = ensemble
    metadata.pop('poses', None)
    run_id = uuid4().hex
    EXPORTS.mkdir(parents=True, exist_ok=True)
    rows, seen = [], {}
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
        for i, record in enumerate(records, 1):
            ext = record['format'] if record['format'] in ('pdb','pdbqt','sdf') else 'txt'
            filename = f'poses/pose_{i:03d}.{ext}'
            archive.writestr(filename, record['text'])
            # Exact coordinate duplicates only; this is not an RMSD clustering claim.
            coords = '\n'.join(line[30:54] for line in record['text'].splitlines()
                               if line.startswith(('ATOM  ', 'HETATM')))
            digest = hashlib.sha256((coords or record['text']).encode()).hexdigest()
            duplicate = seen.get(digest, '')
            seen.setdefault(digest, i)
            rows.append({'sample':i, 'file':filename, 'sampling_attempt':record['attempt'],
                         'score':record.get('energy'), 'score_type':record.get('score_type',''),
                         'exact_duplicate_of':duplicate})
        from .docking_report import score_guide_csv, score_guide_rows, score_meanings_json
        pose_rows = []
        for pose in ensemble:
            normalized_file = pose['file']
            archive.writestr(normalized_file, pose['pdb'])
            engine_name = str(best.get('engine') or '')
            native_details = pose.get('native_details') or {}
            row = {
                'engine': engine_name, 'pose': pose['pose_index'], 'normalized_file': normalized_file,
                'sampling_attempt': pose.get('attempt'),
                'score': pose.get('energy'), 'score_type': pose.get('score_type'),
                'native_score': pose.get('native_score'), 'native_score_type': pose.get('native_score_type'),
                'native_details_json': json.dumps(native_details, sort_keys=True),
                'native_score_meanings_json': score_meanings_json(engine_name, native_details),
                'vina_normalized_score_kcal_mol': pose.get('energy'),
                'vina_score_basis': pose.get('score_type'),
                'confidence': pose.get('confidence'), 'rmsd_to_best_pose': pose.get('rmsd_to_best_pose'),
                'atoms_compared': pose.get('atoms_compared'),
            }
            pose_rows.append(row)
        pose_table = io.StringIO()
        pose_fields = list(pose_rows[0]) if pose_rows else ['pose', 'normalized_file']
        pose_writer = csv.DictWriter(pose_table, fieldnames=pose_fields)
        pose_writer.writeheader(); pose_writer.writerows(pose_rows)
        archive.writestr('poses.csv', pose_table.getvalue())
        # Keep scores.csv as an easy-to-import score table, but include the same
        # complete per-pose data so native metrics cannot be lost on export.
        score_table = io.StringIO()
        score_writer = csv.DictWriter(score_table, fieldnames=pose_fields)
        score_writer.writeheader(); score_writer.writerows(pose_rows)
        archive.writestr('scores.csv', score_table.getvalue())
        guide_rows = score_guide_rows(
            [{'engine': best.get('engine')}],
            [{**pose, 'engine': best.get('engine')} for pose in ensemble],
        )
        archive.writestr('score_guide.csv', score_guide_csv(
            [{'engine': best.get('engine')}],
            [{**pose, 'engine': best.get('engine')} for pose in ensemble],
        ))
        manifest = {'target_id':args[0], 'smiles':args[1], 'engine':best['engine'],
                    'exhaustiveness':args[2] if len(args)>2 else kwargs.get('exhaustiveness',8),
                    'box':best.get('box'), 'pose_count':len(records), 'requested_minimum':requested_poses,
                    'complete':len(records)>=requested_poses, 'warnings':warnings,
                    'best_energy':best['energy'], 'poses':pose_rows,
                    'score_guide': guide_rows}
        archive.writestr('results.json', json.dumps(manifest,indent=2))
        archive.writestr('README.txt',
            f'Real engine outputs; requested at least {requested_poses} poses. Repeated samples can converge to the same pose.\n'
            'Open docking_report.html for the formatted report. poses.csv and scores.csv contain Vina-normalized scores, every returned native score component, and RMSD.\n'
            'Open score_guide.csv for units and plain-language meanings of each score. Check complete/warnings in results.json. No poses are fabricated.\n'
            'DiffDock confidence is not binding energy; its sampled SDF poses are before Vina relaxation.\n'
            'best_pose.pdb is the selected, scored result shown in the app. Files persist until removed from the data directory.\n')
        archive.writestr('best_pose.pdb', best['docked_pdb'])
        from .docking_report import build_docking_report_html
        archive.writestr('docking_report.html', build_docking_report_html({
            'title': f"{best.get('engine', 'MOBIDock')} docking report",
            'target_id': args[0], 'smiles': args[1], 'ligand_label': kwargs.get('ligand_label', 'Docked ligand'),
            'ligand_source': kwargs.get('ligand_source', 'Chosen ligand'),
            'requested_pose_count': requested_poses,
            'exhaustiveness': args[2] if len(args) > 2 else kwargs.get('exhaustiveness', 8),
            'box': best.get('box') or {}, 'props': best.get('props') or {},
            'qed': best.get('qed'), 'synth': best.get('synth') or {},
            'score_note': 'Docking scores are computational estimates, not experimentally measured binding free energies.',
            'engines': [{
                'engine': best.get('engine'), 'display_name': best.get('engine'), 'status': 'completed',
                'best_be_kcal_mol': best.get('energy'),
                'be_source': metadata.get('energy_source'), 'pose_count': len(ensemble),
                'requested_pose_count': requested_poses,
            }],
            'poses': [{**pose, 'engine': best.get('engine'), 'display_name': best.get('engine')}
                      for pose in ensemble],
        }))
    destination = EXPORTS / (run_id+'.zip')
    temporary = destination.with_suffix('.tmp')
    temporary.write_bytes(buffer.getvalue()); temporary.replace(destination)
    best.setdefault('engine_metadata', {})['download'] = {
        'url':f'/api/docking-results/{run_id}', 'pose_count':len(records),
        'complete':len(records)>=requested_poses, 'warnings':warnings}
    return best


@router.post('/api/docking-results/export.xlsx')
async def export_docking_workbook(request: Request):
    """Download an Excel workbook containing every pose currently listed in the app."""
    try:
        payload = await request.json()
    except Exception as exc:
        raise HTTPException(400, 'A JSON payload with a rows array is required.') from exc
    rows = payload.get('rows') if isinstance(payload, dict) else None
    if not isinstance(rows, list):
        raise HTTPException(422, 'The workbook request must include a rows array.')
    if len(rows) > 5000:
        raise HTTPException(413, 'The workbook can contain at most 5,000 pose and run rows.')
    from .docking_report import build_docking_workbook
    try:
        content = build_docking_workbook([row for row in rows if isinstance(row, dict)])
    except Exception as exc:
        raise HTTPException(500, f'Could not create the Excel workbook: {exc}') from exc
    return Response(
        content=content,
        media_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
        headers={'Content-Disposition': 'attachment; filename="mobidock-docking-results.xlsx"'},
    )

@router.get('/api/docking-results/{run_id}')
def download_results(run_id: str):
    if not re.fullmatch(r'[0-9a-f]{32}',run_id):
        raise HTTPException(404,'Result not found')
    path = EXPORTS / (run_id+'.zip')
    if not path.is_file():
        raise HTTPException(404,'Result expired or not found')
    return FileResponse(path, media_type='application/zip', filename=f'docking_{run_id}.zip')
