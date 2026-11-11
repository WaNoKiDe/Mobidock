"""Comparable multi-engine docking runs and downloadable result bundles."""
import csv
import hashlib
import io
import json
import time
import zipfile
from pathlib import Path
from uuid import uuid4

from rdkit import Chem

from . import docking
from .config import DATA_DIR
from .engines import ENGINES
from .engines.poses import first_mol_from_sdf, pdbqt_to_mol, pose_to_pdb
from .engines.rescoring import vina_relax_and_score
from .lipinski import lipinski_filter, rejection_message
from .synth import composite_score, synth_filter, synth_rejection_message
from .results_export import sampling_seed

EXPORTS = DATA_DIR / "docking_exports"
POSE_LIMIT = 20


def _mol_from_pose(text: str, fmt: str):
    if fmt == "sdf":
        return first_mol_from_sdf(text)
    if fmt == "pdbqt":
        return pdbqt_to_mol(text)
    try:
        return Chem.MolFromPDBBlock(text, sanitize=True, removeHs=False, proximityBonding=True)
    except Exception:
        return None


def _raw_pose_records(result: dict, engine_name: str, attempt: int) -> list[dict]:
    metadata = result.get("engine_metadata") or {}
    raw = metadata.get("export_poses") or []
    if raw:
        records = []
        for item in raw:
            text = str(item.get("text") or "")
            if not text.strip():
                continue
            records.append({
                "text": text,
                "format": str(item.get("format") or "pdb"),
                "native_score": item.get("native_score", item.get("energy")),
                "native_score_type": str(item.get("native_score_type") or item.get("score_type") or "engine score"),
                "execution_mode": item.get("execution_mode") or metadata.get("mode") or "",
                "common_be": item.get("common_be"),
                "common_be_source": item.get("common_be_source", ""),
                "common_pose_pdb": item.get("common_pose_pdb", ""),
                "native_details": {
                    **(item.get("native_details") if isinstance(item.get("native_details"), dict) else {}),
                    **{k: item.get(k) for k in ("cnn_score", "cnn_affinity", "karma_score", "karma_score_ff_corrected", "karma_score_aligned") if item.get(k) is not None},
                },
                "attempt": attempt,
            })
        return records

    # KarmaDock currently emits one corrected pose per call. The comparison runner
    # repeats it with independent seeds to match the shared pose request.
    return [{
        "text": result.get("docked_pdbqt") or result.get("docked_pdb") or "",
        "format": "pdbqt" if result.get("docked_pdbqt") else "pdb",
        "native_score": (metadata.get("karma_score_ff_corrected")
                         if metadata.get("karma_score_ff_corrected") is not None
                         else metadata.get("karma_score")),
        "native_score_type": metadata.get("native_score_type") or "KarmaDock score (native units)",
        "execution_mode": metadata.get("mode") or "",
        "native_details": {
            **(metadata.get("native_details") if isinstance(metadata.get("native_details"), dict) else {}),
            **{k: metadata.get(k) for k in ("karma_score", "karma_score_ff_corrected", "karma_score_aligned") if metadata.get(k) is not None},
        },
        "attempt": attempt,
        "common_be": result.get("energy"),
        "common_be_source": metadata.get("energy_source") or "Vina local score after KarmaDock",
    }]


def _decorate_pose(record: dict, engine_name: str, result: dict, receptor_path, center, size) -> dict:
    text, fmt = record["text"], record["format"]
    mol = _mol_from_pose(text, fmt)
    record["pdb"] = record.get("common_pose_pdb") or pose_to_pdb(text)
    record["common_be"] = record.get("common_be")
    record["common_be_source"] = record.get("common_be_source", "")
    if record["common_be"] is None and engine_name in {"vina", "gnina"}:
        # Both engines provide a Vina-family affinity in kcal/mol.
        try:
            record["common_be"] = float(record.get("native_score"))
            record["common_be_source"] = "engine affinity (kcal/mol; Vina-family score)"
        except (TypeError, ValueError):
            pass
    if record["common_be"] is None and engine_name == "diffdock" and mol is not None:
        try:
            energy, relaxed = vina_relax_and_score(Path(receptor_path), mol, center, size)
            record["common_be"] = float(energy)
            record["common_be_source"] = "Vina local relaxation and score (kcal/mol)"
            record["relaxed_pdbqt"] = relaxed
            record["pdb"] = pose_to_pdb(relaxed)
        except Exception as exc:
            record["score_warning"] = str(exc)
    if record["common_be"] is None and engine_name == "karmadock":
        record["common_be"] = result.get("energy")
        record["common_be_source"] = (result.get("engine_metadata") or {}).get("energy_source", "Vina local score (kcal/mol)")
    if not record["pdb"]:
        record["pdb"] = pose_to_pdb(text)
    coords = "\n".join(
        line[30:54] for line in record["pdb"].splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    )
    record["coordinate_hash"] = hashlib.sha256((coords or text).encode("utf-8", errors="replace")).hexdigest()
    return record


def _run_engine(engine_name: str, target_id: str, smiles: str, exhaustiveness: int,
                lipinski_mode: str, synth_mode: str, box: dict, n_poses: int,
                diffdock_api_key: str | None = None) -> dict:
    engine = ENGINES[engine_name]
    started = time.perf_counter()
    problems = engine.problems()
    if problems:
        return {
            "engine": engine_name, "display_name": engine.display_name,
            "status": "not_configured", "message": "; ".join(problems),
            "requested_pose_count": n_poses, "pose_count": 0, "poses": [],
            "runtime_seconds": 0.0,
        }

    pose_records, result_rows, warnings = [], [], []
    # KarmaDock can infer several independently seeded ligand entries in one
    # model load. Other engines return their requested ensemble per invocation.
    for attempt in range(n_poses):
        if len(pose_records) >= n_poses:
            break
        token = sampling_seed.set(42 + attempt * 7919)
        try:
            try:
                result = docking.dock_smiles(
                    target_id, smiles, exhaustiveness, engine_name,
                    lipinski_mode, synth_mode, box_override=box,
                    n_poses=max(1, n_poses - len(pose_records)),
                    diffdock_api_key=diffdock_api_key if engine_name == "diffdock" else None,
                )
            except Exception as exc:
                warnings.append(f"Sampling attempt {attempt + 1} stopped: {type(exc).__name__}: {exc}")
                break
        finally:
            sampling_seed.reset(token)
        result_rows.append(result)
        batch = _raw_pose_records(result, engine_name, attempt + 1)
        if not batch:
            batch = _raw_pose_records({"docked_pdbqt": result.get("docked_pdbqt"),
                                       "docked_pdb": result.get("docked_pdb"),
                                       "engine_metadata": result.get("engine_metadata")}, engine_name, attempt + 1)
        for record in batch:
            if len(pose_records) >= n_poses:
                break
            receptor_path = docking.RECEPTOR_PATHS.get(target_id)
            pose_records.append(_decorate_pose(record, engine_name, result, receptor_path,
                                                tuple(box["center"]), tuple(box["size"])))
        if len(pose_records) >= n_poses:
            break

    if not pose_records:
        raise RuntimeError("The engine returned no usable poses.")
    valid = [p for p in pose_records if isinstance(p.get("common_be"), (int, float))]
    best_pose = min(valid, key=lambda p: float(p["common_be"])) if valid else pose_records[0]
    best_result = min(result_rows, key=lambda r: float(r.get("energy", 1e9))) if result_rows else {}
    for rank, pose in enumerate(sorted(valid, key=lambda p: float(p["common_be"])), 1):
        pose["be_rank"] = rank
    for rank, pose in enumerate(pose_records, 1):
        pose["engine_rank"] = rank

    seen = {}
    for pose in pose_records:
        pose["exact_duplicate_of"] = seen.get(pose["coordinate_hash"])
        seen.setdefault(pose["coordinate_hash"], pose["engine_rank"])
    raw_meta = best_result.get("engine_metadata") or {}
    meta = {k: v for k, v in raw_meta.items() if k not in {"export_poses", "poses"}}
    return {
        "engine": engine_name, "display_name": engine.display_name,
        "status": "completed", "message": "; ".join(warnings),
        "requested_pose_count": n_poses, "pose_count": len(pose_records),
        "pose_count_complete": len(pose_records) == n_poses,
        "best_be_kcal_mol": best_pose.get("common_be"),
        "be_source": best_pose.get("common_be_source") or "unavailable",
        "best_native_score": best_pose.get("native_score"),
        "native_score_type": best_pose.get("native_score_type", "engine score"),
        "native_details": best_pose.get("native_details", {}),
        "record_energy": best_result.get("energy"),
        "mode": meta.get("mode", ""), "engine_metadata": meta,
        "runtime_seconds": round(time.perf_counter() - started, 2),
        "best_pose_pdb": best_pose.get("pdb", ""),
        "props": best_result.get("props", {}), "qed": best_result.get("qed"),
        "synth": best_result.get("synth", {}),
        "poses": pose_records,
    }


def run_comparison(target_id: str, smiles: str, exhaustiveness: int = 8,
                   lipinski_mode: str = "standard", synth_mode: str = "standard",
                   box_override: dict | None = None, n_poses: int = 10,
                   ligand_source: str = "chosen", ligand_label: str = "Chosen ligand",
                   player_name: str = "Anonymous", engine_name: str | None = None,
                   diffdock_api_key: str | None = None) -> dict:
    """Dock with every registered engine or one selected engine and build a portable ZIP report."""
    if target_id not in docking.BOX_PARAMS:
        raise ValueError(f"Unknown target: {target_id!r}.")
    smiles = (smiles or "").strip()
    base_mol = Chem.MolFromSmiles(smiles) if smiles else None
    if base_mol is None:
        raise ValueError("Enter a valid ligand SMILES before docking.")
    lipinski = lipinski_filter(base_mol, mode=lipinski_mode)
    if not lipinski["passed"]:
        raise ValueError(rejection_message(lipinski))
    synth_mode = (synth_mode or "standard").strip().lower()
    synth = synth_filter(base_mol, mode=synth_mode)
    if not synth["passed"]:
        raise ValueError(synth_rejection_message(synth))
    n_poses = max(1, min(POSE_LIMIT, int(n_poses or 10)))
    box = box_override or docking.BOX_PARAMS[target_id]
    requested_engines = [str(engine_name or "").strip().lower()] if engine_name else list(ENGINES)
    unknown = [name for name in requested_engines if name not in ENGINES]
    if unknown:
        raise ValueError(f"Unknown docking engine: {unknown[0]!r}. Choose from {list(ENGINES)}")
    player_name = (player_name or "").strip() or "Anonymous"
    summary, engine_results, all_poses = [], [], []
    for engine in requested_engines:
        try:
            row = _run_engine(engine, target_id, smiles, exhaustiveness,
                              lipinski_mode, synth_mode, box, n_poses,
                              diffdock_api_key if engine == "diffdock" else None)
        except Exception as exc:
            row = {
                "engine": engine, "display_name": ENGINES[engine].display_name,
                "status": "failed", "message": f"{type(exc).__name__}: {exc}",
                "requested_pose_count": n_poses, "pose_count": 0,
                "pose_count_complete": False, "runtime_seconds": 0.0, "poses": [],
            }
        for pose in row.get("poses", []):
            pose["engine"] = engine
            all_poses.append(pose)
        if row.get("status") == "completed" and row.get("record_energy") is not None:
            try:
                from . import recorddb
                synth_result = row.get("synth") or {}
                synth_norm = synth_result.get("norm") if synth_result.get("available") else None
                score = composite_score(float(row["record_energy"]), float(row.get("qed") or 0), synth_norm)
                row["history_saved"] = bool(recorddb.save_docking(
                    target_id, engine, exhaustiveness, player_name, smiles,
                    row["record_energy"], row.get("qed"), score,
                    row.get("best_pose_pdb") or "", row.get("props") or {},
                    synth_result, box))
            except Exception as exc:
                row["history_saved"] = False
                message = str(row.get("message") or "").strip()
                row["message"] = (message + "; " if message else "") + f"Docking-history save failed: {exc}"
        else:
            row["history_saved"] = False
        summary.append({k: v for k, v in row.items() if k not in {"poses", "best_pose_pdb", "props", "synth", "record_energy"}})
        engine_results.append({k: v for k, v in row.items() if k not in {"poses", "record_energy"}})

    # Calculate each engine's pose RMSDs against that engine's lowest common
    # binding-energy pose in the original receptor coordinate frame.
    try:
        from .poses import pose_mol, rmsd_common
        by_engine = {}
        for pose in all_poses:
            by_engine.setdefault(pose.get("engine"), []).append(pose)
        for engine_poses in by_engine.values():
            scored = [pose for pose in engine_poses if isinstance(pose.get("common_be"), (int, float))]
            if not scored:
                continue
            reference = min(scored, key=lambda pose: float(pose["common_be"]))
            reference_mol = pose_mol(reference.get("pdb") or "")
            reference["rmsd_to_best_pose"] = 0.0
            reference["atoms_compared"] = None
            for pose in engine_poses:
                if pose is reference:
                    continue
                compared = rmsd_common(reference_mol, pose_mol(pose.get("pdb") or ""), min_atoms=1)
                pose["rmsd_to_best_pose"] = round(float(compared[0]), 3) if compared else None
                pose["atoms_compared"] = int(compared[1]) if compared else None
    except Exception:
        pass

    run_id = uuid4().hex
    EXPORTS.mkdir(parents=True, exist_ok=True)
    archive_path = EXPORTS / f"{run_id}.zip"
    from .docking_report import score_guide_csv, score_guide_rows, score_meanings_json
    for pose in all_poses:
        details = pose.get("native_details") or {}
        pose["vina_normalized_score_kcal_mol"] = pose.get("common_be")
        pose["vina_score_basis"] = pose.get("common_be_source")
        pose["native_score_meanings_json"] = score_meanings_json(pose.get("engine"), details)
    summary_fields = ["engine", "display_name", "status", "message", "requested_pose_count", "pose_count",
                      "pose_count_complete", "best_be_kcal_mol", "be_source", "best_native_score",
                      "native_score_type", "best_native_details_json", "best_native_score_meanings_json",
                      "history_saved", "mode", "runtime_seconds"]
    pose_fields = ["engine", "engine_rank", "be_rank", "common_be", "common_be_source", "native_score",
                   "vina_normalized_score_kcal_mol", "vina_score_basis", "native_score_type", "native_details", "native_details_json",
                   "native_score_meanings_json", "attempt", "format", "exact_duplicate_of",
                   "execution_mode", "score_warning", "rmsd_to_best_pose", "atoms_compared", "file", "normalized_file"]
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED) as archive:
        out = io.StringIO(newline="")
        writer = csv.DictWriter(out, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        for item in summary:
            summary_row = {key: value for key, value in item.items() if key in summary_fields}
            summary_row["best_native_details_json"] = json.dumps(item.get("native_details") or {}, sort_keys=True)
            summary_row["best_native_score_meanings_json"] = score_meanings_json(item.get("engine"), item.get("native_details") or {})
            writer.writerow(summary_row)
        archive.writestr("engine_summary.csv", out.getvalue())
        pose_rows = []
        for pose in all_poses:
            engine = pose.get("engine", "engine")
            ext = pose.get("format") if pose.get("format") in {"pdb", "pdbqt", "sdf"} else "txt"
            name = f"poses/{engine}/pose_{int(pose.get('engine_rank', 0)):03d}.{ext}"
            archive.writestr(name, pose["text"])
            pdb_text = pose.get("pdb") or pose_to_pdb(pose["text"])
            pose["file"] = name
            if pdb_text:
                normalized_name = f"poses/{engine}/pose_{int(pose.get('engine_rank', 0)):03d}_normalized.pdb"
                archive.writestr(normalized_name, pdb_text)
                pose["normalized_file"] = normalized_name
            else:
                pose["normalized_file"] = ""
            native_details = pose.get("native_details") or {}
            row = {k: pose.get(k) for k in pose_fields if k != "file"}
            row["vina_normalized_score_kcal_mol"] = pose.get("common_be")
            row["vina_score_basis"] = pose.get("common_be_source")
            row["native_details"] = json.dumps(native_details, sort_keys=True)
            row["native_details_json"] = json.dumps(native_details, sort_keys=True)
            row["native_score_meanings_json"] = score_meanings_json(pose.get("engine"), native_details)
            row["file"] = name
            pose_rows.append(row)
        out = io.StringIO(newline="")
        writer = csv.DictWriter(out, fieldnames=pose_fields, extrasaction="ignore")
        writer.writeheader()
        for row in pose_rows:
            writer.writerow(row)
        archive.writestr("poses.csv", out.getvalue())
        for result in engine_results:
            if result.get("best_pose_pdb"):
                archive.writestr(f"best_pose_{result['engine']}.pdb", result["best_pose_pdb"])
        guide_rows = score_guide_rows(summary, all_poses)
        manifest = {
            "target_id": target_id, "ligand_source": ligand_source, "ligand_label": ligand_label,
            "smiles": smiles, "exhaustiveness": exhaustiveness,
            "requested_pose_count_per_engine": n_poses,
            "box": {"center": list(box["center"]), "size": list(box["size"])},
            "score_note": "BE values are docking scores in kcal/mol, not experimentally measured binding free energies. Vina-family affinities are compared as an approximate common scale; native engine scores remain separate.",
            "engines": summary,
            "poses": pose_rows,
            "score_guide": guide_rows,
        }
        archive.writestr("results.json", json.dumps(manifest, indent=2, default=str))
        archive.writestr("score_guide.csv", score_guide_csv(summary, all_poses))
        from .docking_report import build_docking_report_html
        profile = next((result for result in engine_results if result.get("props")), {})
        report_poses = []
        for pose in all_poses:
            report_poses.append({
                **pose,
                "display_name": ENGINES.get(pose.get("engine"), ENGINES[requested_engines[0]]).display_name,
                "pose_index": pose.get("engine_rank"), "energy": pose.get("common_be"),
                "score_type": pose.get("common_be_source"), "confidence": pose.get("confidence"),
            })
        archive.writestr("docking_report.html", build_docking_report_html({
            "title": "MOBIDock docking report",
            "target_id": target_id, "smiles": smiles, "ligand_label": ligand_label,
            "ligand_source": ligand_source, "requested_pose_count": n_poses,
            "exhaustiveness": exhaustiveness,
            "box": {"center": list(box["center"]), "size": list(box["size"])},
            "props": profile.get("props") or {}, "qed": profile.get("qed"),
            "synth": profile.get("synth") or synth, "score_note": manifest["score_note"],
            "engines": summary, "poses": report_poses,
        }))
        archive.writestr("README.txt",
            "MOBIDock docking results\n\n"
            "Open docking_report.html for the formatted report. engine_summary.csv compares returned poses and the best reported kcal/mol docking score. "
            "poses.csv reports every pose, the Vina-normalized score, all returned engine-native score components and meanings, RMSD, and coordinate file.\n"
            "engine_summary.csv includes each engine's best normalized and native scores. score_guide.csv explains score units and meanings.\n"
            "Pose counts are real outputs and are never padded with copies. If an engine returns fewer "
            "than requested, inspect pose_count_complete and the engine message.\n"
            "Every completed engine run is also saved to MOBIDock docking history.\n"
            "Docking scores are predictions, not experimentally measured binding free energies. "
            "Vina-family scores are approximate cross-engine comparisons, not guaranteed rankings.\n")

    return {
        "run_id": run_id, "target_id": target_id, "smiles": smiles,
        "synth": synth,
        "ligand_source": ligand_source, "ligand_label": ligand_label,
        "requested_pose_count": n_poses, "requested_engines": requested_engines,
        "run_mode": "single" if engine_name else "all",
        "score_note": "BE values are docking scores in kcal/mol, not experimental binding free energies. Vina-family scores are an approximate common scale; each engine's native score is shown separately.",
        "summary": summary, "results": engine_results,
        "pose_table": [{k: p.get(k) for k in pose_fields if k != "file"} for p in all_poses],
        "download": {"url": f"/api/docking-results/{run_id}",
                     "filename": f"mobidock_{requested_engines[0]}_{run_id[:8]}.zip" if engine_name else f"mobidock_comparison_{run_id[:8]}.zip"},
    }
