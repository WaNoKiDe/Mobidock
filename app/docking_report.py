"""Self-contained, printable HTML reports included with docking result bundles."""
from datetime import datetime, timezone
from html import escape
import json
import re


def _text(value, fallback="—"):
    if value is None or value == "":
        return fallback
    return escape(str(value))


def _number(value, digits=2, suffix=""):
    try:
        number = float(value)
        if number != number or number in (float("inf"), float("-inf")):
            return "—"
        return f"{number:.{digits}f}{suffix}"
    except (TypeError, ValueError):
        return "—"


_SCORE_DEFINITIONS = {
    "minimizedaffinity": ("kcal/mol", "GNINA locally relaxes the pose, then estimates pocket fit using a Vina-family score with contact rewards and clash/strain penalties. Lower (more negative) is more favorable; this is a model estimate, not measured binding energy."),
    "vinabindingaffinity": ("kcal/mol", "AutoDock Vina combines steric, hydrophobic, hydrogen-bond and flexibility terms to estimate how well the ligand fits the pocket. Lower (more negative) is more favorable; it is not an experimental free energy."),
    "cnnscore": ("0–1", "GNINA's neural network estimates whether a pose is close to a correct, near-native pose (roughly within 2 Å). Higher means the pose looks more plausible; it does not mean stronger binding."),
    "cnnaffinity": ("pK", "GNINA predicts affinity on a logarithmic pK scale related to how tightly a ligand binds. Higher predicts tighter binding; this is a model prediction, not an assay result."),
    "cnnvs": ("CNNaffinity × CNNscore", "GNINA's virtual-screening rank combines predicted affinity with pose plausibility. Higher ranks more favorably; it is a composite score, not a physical energy."),
    "cnnaffinityvariance": ("pK²", "Spread among GNINA's affinity predictions. Lower means the predictions agree more closely; this describes uncertainty, not binding strength."),
    "karmascore": ("KarmaDock native units", "KarmaDock's mixture-density model assigns support to the pose and its predicted interaction strength. Higher means the model favors it more; this is not kcal/mol."),
    "karmascoreff": ("KarmaDock native units", "The KarmaDock model score after force-field pose correction, reflecting how much the model supports the corrected interaction. Higher is more favorable within KarmaDock; it is not kcal/mol."),
    "karmascoreffcorrected": ("KarmaDock native units", "The KarmaDock model score after force-field pose correction, reflecting how much the model supports the corrected interaction. Higher is more favorable within KarmaDock; it is not kcal/mol."),
    "karmascorealigned": ("KarmaDock native units", "The KarmaDock model score after pose-alignment correction. Higher means the model supports that aligned pose more; it is not kcal/mol."),
    "confidence": ("model score (unitless)", "DiffDock's confidence model estimates how plausible the predicted pose geometry is. Higher means more confidence in the pose, not stronger binding or lower energy."),
    "nativescore": ("engine-native units", "The engine's own primary measure: GNINA and Vina use an energy-like fit score (lower is better), DiffDock reports pose confidence (higher is better), and KarmaDock reports model support (higher is better). These scales are not comparable across engines."),
    "vinanormalizedscorekcalmol": ("kcal/mol", "Vina locally relaxes and scores each pose using pocket contacts and clash/strain penalties, providing a shared comparison scale. Lower (more negative) ranks as more favorable; it is not experimental binding energy."),
    "commonbe": ("kcal/mol", "Vina locally relaxes and scores each pose using pocket contacts and clash/strain penalties, providing a shared comparison scale. Lower (more negative) ranks as more favorable; it is not experimental binding energy."),
}


def _score_key(name):
    return re.sub(r"[^a-z0-9]", "", str(name or "").lower())


def score_definition(engine, field):
    """Return a portable unit and plain-language explanation for score fields."""
    key = _score_key(field)
    if key in _SCORE_DEFINITIONS:
        return _SCORE_DEFINITIONS[key]
    engine_name = str(engine or "docking engine")
    engine_key = str(engine or "").strip().lower()
    if key in {"native", "nativescore"}:
        if engine_key == "gnina":
            return ("kcal/mol", "GNINA's primary minimizedAffinity is a Vina-family pocket-fit score after local pose relaxation. Lower (more negative) is more favorable; it is not experimental binding energy.")
        if engine_key == "vina":
            return _SCORE_DEFINITIONS["vinabindingaffinity"]
        if engine_key == "diffdock":
            return _SCORE_DEFINITIONS["confidence"]
        if engine_key == "karmadock":
            return _SCORE_DEFINITIONS["karmascore"]
        return ("engine-native units", f"The primary score reported by {engine_name}; its scale and favorable direction are engine-specific, so do not compare it directly across engines.")
    if "score" in key or "affinity" in key or "energy" in key:
        if engine_key in {"gnina", "vina"}:
            return ("as reported by engine", f"This {engine_name} component is an energy-like fit score; lower is usually more favorable. Its precise terms and units depend on the field name, and it is not an experimental measurement.")
        if engine_key in {"diffdock", "karmadock"}:
            return ("as reported by engine", f"This {engine_name} component is a model score, not physical energy; higher generally means the model ranks the pose more favorably. Keep it in its native scale.")
        return ("as reported by engine", f"This component stays in {engine_name}'s native scale; the favorable direction depends on the particular metric, so it is not directly comparable across engines.")
    return ("as reported by engine", f"Additional output from {engine_name}; its direction depends on the metric and it is not a binding-energy measurement.")


def score_guide_rows(engines, poses):
    """List the normalized score and every native score field present in a run."""
    engine_names = {str(row.get("engine") or "").lower() for row in (engines or []) if isinstance(row, dict)}
    detail_keys = {}
    for pose in poses or []:
        if not isinstance(pose, dict):
            continue
        engine = str(pose.get("engine") or "").lower()
        details = pose.get("native_details") if isinstance(pose.get("native_details"), dict) else {}
        detail_keys.setdefault(engine, set()).update(str(key) for key in details)
    rows = [{
        "engine": "All engines", "score_field": "vina_normalized_score_kcal_mol",
        "units": _SCORE_DEFINITIONS["vinanormalizedscorekcalmol"][0],
        "meaning": _SCORE_DEFINITIONS["vinanormalizedscorekcalmol"][1],
    }]
    for engine in sorted(engine_names | set(detail_keys)):
        names = set(detail_keys.get(engine, set()))
        names.add("native_score")
        if engine == "gnina":
            names.update(("minimizedAffinity", "CNNscore", "CNNaffinity", "CNN_VS", "CNNaffinity_variance"))
        elif engine == "karmadock":
            names.update(("karma_score", "karma_score_ff", "karma_score_aligned"))
        elif engine == "diffdock":
            names.add("confidence")
        elif engine == "vina":
            names.add("Vina binding affinity")
        for name in sorted(names, key=str.lower):
            units, meaning = score_definition(engine, name)
            rows.append({"engine": engine or "Engine", "score_field": name, "units": units, "meaning": meaning})
    return rows


def score_meanings_json(engine, details):
    details = details if isinstance(details, dict) else {}
    fields = {"native_score": score_definition(engine, "native_score")[0:2]}
    engine_key = str(engine or "").lower()
    if engine_key == "gnina":
        fields["minimizedAffinity"] = score_definition(engine, "minimizedAffinity")
    elif engine_key == "diffdock":
        fields["confidence"] = score_definition(engine, "confidence")
    elif engine_key == "vina":
        fields["Vina binding affinity"] = score_definition(engine, "Vina binding affinity")
    elif engine_key == "karmadock":
        fields["karma_score_ff"] = score_definition(engine, "karma_score_ff")
    fields.update({
        str(name): score_definition(engine, name)
        for name in sorted(details, key=str.lower)
    })
    return json.dumps({
        str(name): {"units": definition[0], "meaning": definition[1]}
        for name, definition in fields.items()
    }, sort_keys=True)


def _native_details_html(engine, details):
    details = details if isinstance(details, dict) else {}
    if not details:
        return ""
    items = []
    for name, value in sorted(details.items(), key=lambda item: str(item[0]).lower()):
        units, meaning = score_definition(engine, name)
        items.append(
            f"<li><b>{escape(str(name))}</b>: {_text(value)} <span class='unit'>({escape(units)})</span>"
            f"<br><span class='foot'>{escape(meaning)}</span></li>"
        )
    return "<ul class='native-score-list'>" + "".join(items) + "</ul>"


def score_guide_csv(engines, poses):
    """Return an RFC-4180 CSV text listing score names, units and meanings."""
    import csv
    import io
    output = io.StringIO(newline="")
    rows = score_guide_rows(engines, poses)
    writer = csv.DictWriter(output, fieldnames=["engine", "score_field", "units", "meaning"])
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _molecule_svg(smiles):
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = Chem.MolFromSmiles(str(smiles or ""))
        if mol is None:
            return "<div class='structure-empty'>Structure drawing unavailable</div>"
        drawer = rdMolDraw2D.MolDraw2DSVG(520, 260)
        drawer.drawOptions().clearBackground = False
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return "<div class='structure-empty'>Structure drawing unavailable</div>"


def _property_rows(props):
    props = props if isinstance(props, dict) else {}
    fields = [
        ("Molecular weight", props.get("mw"), " Da"),
        ("LogP", props.get("logp"), ""),
        ("Hydrogen-bond donors", props.get("hbd"), ""),
        ("Hydrogen-bond acceptors", props.get("hba"), ""),
        ("Topological polar surface area", props.get("tpsa"), " Å²"),
        ("Rotatable bonds", props.get("rot_bonds"), ""),
    ]
    rows = "".join(
        f"<tr><th>{escape(label)}</th><td>{_number(value, 2)}{escape(unit)}</td></tr>"
        for label, value, unit in fields
    )
    lipinski = props.get("lipinski1") if isinstance(props.get("lipinski1"), dict) else {}
    lipinski_text = "Lipinski screen not recorded"
    if lipinski:
        violations = lipinski.get("violations")
        if isinstance(violations, list):
            violations = len(violations)
        lipinski_text = (
            f"{_text(lipinski.get('passed'), 'Not recorded')} · "
            f"{_text(violations, '—')} violation(s) · mode {_text(lipinski.get('mode'))}"
        )
    rows += f"<tr><th>Lipinski screen</th><td>{lipinski_text}</td></tr>"
    return rows


def _engine_rows(engines):
    rows = []
    for engine in engines or []:
        rows.append(
            "<tr>"
            f"<td>{_text(engine.get('display_name') or engine.get('engine'))}</td>"
            f"<td>{_text(engine.get('status'))}</td>"
            f"<td>{_number(engine.get('best_be_kcal_mol', engine.get('energy')), 3, ' kcal/mol')}</td>"
            f"<td>{_text(engine.get('be_source') or engine.get('energy_source') or engine.get('native_score_type'))}</td>"
            f"<td>{_text(engine.get('pose_count'))} / {_text(engine.get('requested_pose_count'))}</td>"
            f"<td>{_number(engine.get('runtime_seconds'), 2, ' s')}</td>"
            "</tr>"
        )
    return "".join(rows) or "<tr><td colspan='6'>No engine summary was returned.</td></tr>"


def _pose_rows(poses):
    rows = []
    for pose in poses or []:
        if not isinstance(pose, dict):
            continue
        score = pose.get("common_be", pose.get("energy"))
        native = pose.get("native_score")
        rmsd = pose.get("rmsd_to_best_pose", pose.get("rmsd"))
        rows.append(
            "<tr>"
            f"<td>{_text(pose.get('display_name') or pose.get('engine'))}</td>"
            f"<td>{_text(pose.get('pose_index') or pose.get('engine_rank'))}</td>"
            f"<td>{_number(score, 3)}</td>"
            f"<td>{_text(pose.get('common_be_source') or pose.get('score_type') or pose.get('energy_source'))}</td>"
            f"<td>{_number(native, 3)}</td>"
            f"<td>{_text(pose.get('native_score_type'))}{_native_details_html(pose.get('engine'), pose.get('native_details'))}</td>"
            f"<td>{_number(rmsd, 3, ' Å')}</td>"
            f"<td>{_text(pose.get('atoms_compared'))}</td>"
            f"<td>{_number(pose.get('confidence'), 3)}</td>"
            f"<td class='file'>{_text(pose.get('normalized_file') or pose.get('file'))}</td>"
            "</tr>"
        )
    return "".join(rows) or "<tr><td colspan='10'>No pose-level results were returned.</td></tr>"


def _fingerprint_section(smiles):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, BRICS, MACCSkeys

        mol = Chem.MolFromSmiles(str(smiles or ""))
        if mol is None:
            raise ValueError("invalid SMILES")
        bit_info = {}
        fingerprint = AllChem.GetMorganFingerprint(mol, 2, bitInfo=bit_info)
        counts = fingerprint.GetNonzeroElements()
        morgan_rows = []
        for bit, count in sorted(counts.items(), key=lambda item: item[0])[:100]:
            atom_index, radius = bit_info[bit][0]
            bonds = list(Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index))
            atom_ids = {atom_index}
            for bond_index in bonds:
                bond = mol.GetBondWithIdx(bond_index)
                atom_ids.update((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            fragment = Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(atom_ids), bondsToUse=bonds,
                                                canonical=True, isomericSmiles=True)
            morgan_rows.append(
                f"<tr><td>{bit}</td><td>{count}</td><td>{radius}</td><td class='file'>{escape(fragment)}</td></tr>"
            )
        maccs = list(MACCSkeys.GenMACCSKeys(mol).GetOnBits())
        brics = sorted(BRICS.BRICSDecompose(mol))
        morgan_table = "".join(morgan_rows) or "<tr><td colspan='4'>No Morgan bits were returned.</td></tr>"
        brics_text = " · ".join(escape(fragment) for fragment in brics) or "No BRICS fragments identified."
        return f"""<div class='metric fingerprint-metrics'><div><span>Morgan radius</span><b>2</b></div>
<div><span>Unique Morgan bits</span><b>{len(counts)}</b></div><div><span>MACCS keys set</span><b>{len(maccs)}</b></div></div>
<p><b>MACCS key indices:</b> <span class='file'>{', '.join(str(bit) for bit in maccs) or 'None'}</span></p>
<div class='card'><table><thead><tr><th>Morgan bit</th><th>Count</th><th>Radius</th><th>Representative fragment SMILES</th></tr></thead><tbody>{morgan_table}</tbody></table>
<p class='foot'>Showing up to 100 unique Morgan bits. A bit can represent more than one atom environment.</p></div>
<h3>BRICS fragments</h3><div class='smiles'>{brics_text}</div>"""
    except Exception:
        return "<div class='card'>Fingerprint details are unavailable for this ligand.</div>"


def build_docking_report_html(report):
    """Return a self-contained report suitable for inclusion in a ZIP archive."""
    report = report if isinstance(report, dict) else {}
    props = report.get("props") if isinstance(report.get("props"), dict) else {}
    synth = report.get("synth") if isinstance(report.get("synth"), dict) else {}
    engines = report.get("engines") or []
    poses = report.get("poses") or []
    fingerprints = _fingerprint_section(report.get("smiles"))
    box = report.get("box") if isinstance(report.get("box"), dict) else {}
    center, size = box.get("center"), box.get("size")
    box_text = "—"
    if isinstance(center, (list, tuple)) and isinstance(size, (list, tuple)):
        box_text = (
            "Center (Å): " + ", ".join(_number(v, 2) for v in center)
            + " · Size (Å): " + " × ".join(_number(v, 1) for v in size)
        )
    qed = _number(report.get("qed"), 3)
    sa_text = "Not available"
    if synth.get("available"):
        sa_text = (
            f"SA {_number(synth.get('sa'), 2)} / 10 · {_text(synth.get('tier_label'))} · "
            f"{_number(synth.get('score'), 0)} / 100"
        )
    elif synth:
        sa_text = _text(synth.get("error"), "Scorer unavailable")
    note = _text(
        report.get("score_note"),
        "Docking scores are computational estimates, not experimentally measured binding free energies.",
    )
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    title = _text(report.get("title"), "MOBIDock docking report")
    smiles = _text(report.get("smiles"), "—")
    engine_count = len(engines)
    pose_count = len(poses)
    score_guide = score_guide_rows(engines, poses)
    score_guide_html = "".join(
        f"<tr><td>{_text(row['engine'])}</td><td class='file'>{_text(row['score_field'])}</td>"
        f"<td>{_text(row['units'])}</td><td>{_text(row['meaning'])}</td></tr>"
        for row in score_guide
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>
:root{{color-scheme:light;--ink:#172033;--muted:#526078;--line:#d8dfeb;--navy:#15345b;--blue:#eaf2fc;--teal:#dff5f0;--paper:#fff}}
*{{box-sizing:border-box}}body{{margin:0;background:#eef2f7;color:var(--ink);font:14px/1.45 Arial,Helvetica,sans-serif}}
.page{{max-width:1040px;margin:28px auto;padding:0 18px}}.sheet{{background:var(--paper);padding:34px 42px;border:1px solid var(--line);border-radius:10px;box-shadow:0 6px 24px #10233b12}}
.masthead{{display:flex;justify-content:space-between;gap:20px;border-bottom:3px solid var(--navy);padding-bottom:16px;margin-bottom:22px}}
h1{{font-size:25px;line-height:1.15;margin:0;color:var(--navy)}}.kicker{{color:var(--muted);font-size:11px;font-weight:700;letter-spacing:.13em;text-transform:uppercase;margin-bottom:7px}}.date{{color:var(--muted);white-space:nowrap;font-size:12px}}
h2{{font-size:16px;color:var(--navy);margin:25px 0 10px;border-bottom:1px solid var(--line);padding-bottom:6px}}
.grid{{display:grid;grid-template-columns:1.1fr .9fr;gap:18px}}.card{{border:1px solid var(--line);border-radius:8px;padding:14px 16px;background:#fff}}
.structure{{display:flex;align-items:center;justify-content:center;min-height:220px;background:#fbfcfe;border-radius:6px}}.structure svg{{max-width:100%;height:auto}}
.meta{{display:grid;grid-template-columns:130px 1fr;gap:7px 12px;margin:0}}.meta dt{{color:var(--muted);font-weight:700}}.meta dd{{margin:0;overflow-wrap:anywhere}}
.smiles{{font:12px/1.45 Consolas,monospace;background:#f5f7fa;padding:10px;border-radius:5px;overflow-wrap:anywhere}}
table{{width:100%;border-collapse:collapse;font-size:12px}}th,td{{padding:8px 9px;border:1px solid var(--line);text-align:left;vertical-align:top}}thead th{{background:var(--navy);color:#fff;font-weight:700}}tbody tr:nth-child(even){{background:#f7f9fc}}
.properties th{{width:48%;background:#f4f7fb;color:var(--muted)}}.metric{{display:grid;grid-template-columns:repeat(3,1fr);gap:9px}}.metric div{{padding:12px;background:var(--blue);border-radius:7px}}.metric span{{display:block;color:var(--muted);font-size:11px}}.metric b{{display:block;color:var(--navy);font-size:16px;margin-top:4px}}
.fingerprint-metrics{{margin:12px 0}}h3{{font-size:13px;color:var(--navy);margin:16px 0 8px}}
.notice{{border-left:4px solid #e7a92f;background:#fff8e6;padding:11px 13px;color:#594b2c;border-radius:4px}}.foot{{margin-top:24px;color:var(--muted);font-size:11px;border-top:1px solid var(--line);padding-top:12px}}.file{{font:11px Consolas,monospace;overflow-wrap:anywhere}}.native-score-list{{padding-left:16px;margin:6px 0 0}}.native-score-list li{{margin:4px 0}}.unit{{color:var(--muted);font-size:10px}}
@media(max-width:720px){{.sheet{{padding:22px 18px}}.grid{{grid-template-columns:1fr}}.masthead{{display:block}}.date{{margin-top:10px}}}}
@media print{{body{{background:#fff}}.page{{margin:0;max-width:none;padding:0}}.sheet{{border:0;border-radius:0;box-shadow:none;padding:14mm}}h2,.card,.notice{{break-inside:avoid}}table{{font-size:9px}}th,td{{padding:5px}}thead{{display:table-header-group}}}}
</style></head><body><main class="page"><article class="sheet">
<header class="masthead"><div><div class="kicker">MOBIDock · Computational docking report</div><h1>{title}</h1></div><div class="date">Generated {generated}</div></header>
<section class="grid"><div class="card"><div class="kicker">Ligand structure</div><div class="structure">{_molecule_svg(report.get('smiles'))}</div><div class="smiles">{smiles}</div></div>
<div class="card"><div class="kicker">Run details</div><dl class="meta">
<dt>Target</dt><dd>{_text(report.get('target_id'))}</dd><dt>Ligand</dt><dd>{_text(report.get('ligand_label'), smiles)}</dd>
<dt>Ligand source</dt><dd>{_text(report.get('ligand_source'))}</dd><dt>Docking engine(s)</dt><dd>{engine_count}</dd>
<dt>Requested poses / engine</dt><dd>{_text(report.get('requested_pose_count'))}</dd><dt>Pose records</dt><dd>{pose_count}</dd>
<dt>Grid box</dt><dd>{escape(box_text)}</dd><dt>Exhaustiveness</dt><dd>{_text(report.get('exhaustiveness'))}</dd>
</dl></div></section>
<section><h2>1 · Molecule profile</h2><div class="grid"><div class="card"><table class="properties"><tbody>{_property_rows(props)}</tbody></table></div>
<div class="metric"><div><span>QED</span><b>{qed}</b></div><div><span>Synthesizability</span><b>{sa_text}</b></div><div><span>Engine runs</span><b>{engine_count}</b></div></div></div></section>
<section><h2>2 · Docking summary</h2><div class="card"><table><thead><tr><th>Engine</th><th>Status</th><th>Best common score</th><th>Score basis</th><th>Returned / requested</th><th>Runtime</th></tr></thead><tbody>{_engine_rows(engines)}</tbody></table></div></section>
<section><h2>3 · Fingerprints and fragments</h2>{fingerprints}</section>
<section><h2>4 · Pose-level results</h2><div class="card"><table><thead><tr><th>Engine</th><th>Pose</th><th>Vina-normalized score (kcal/mol)</th><th>Score basis</th><th>Native score</th><th>Native score type and components</th><th>RMSD to best pose</th><th>Atoms</th><th>Confidence</th><th>Pose file</th></tr></thead><tbody>{_pose_rows(poses)}</tbody></table></div></section>
<section><h2>5 · Score guide</h2><div class="card"><table><thead><tr><th>Engine</th><th>Score field</th><th>Units</th><th>Physical or intuitive meaning · favorable direction</th></tr></thead><tbody>{score_guide_html}</tbody></table></div><p class="foot">GNINA normally returns minimizedAffinity, CNNscore, CNNaffinity and CNN_VS. When a metric is absent from a pose file it is left blank; it is never estimated or fabricated. DiffDock confidence and KarmaDock native scores are separate from the Vina-normalized kcal/mol value.</p></section>
<section><h2>Interpretation and files</h2><div class="notice">{note} RMSD is calculated in the protein coordinate frame against the lowest-energy returned pose when an energy is available. Scores from different engines are approximate and should not be treated as experimental affinity measurements.</div></section>
<footer class="foot">This report summarizes the docking run and its returned poses. The ZIP archive also includes the original and normalized pose-coordinate files, machine-readable CSV tables, and a JSON manifest.</footer>
</article></main></body></html>"""


def build_docking_workbook(pose_rows):
    """Create a dependency-free Excel .xlsx workbook from pose rows."""
    from .spreadsheet_export import build_workbook

    poses = list(pose_rows or [])
    pose_headers = [
        ("Target", "target_id"), ("Ligand", "ligand_label"), ("SMILES", "smiles"),
        ("Engine", "display_name"), ("Engine key", "engine"), ("Run", "batch_id"),
        ("Status", "status"), ("Pose", "pose_index"), ("Returned poses", "pose_count"),
        ("Requested poses", "requested_pose_count"), ("Vina-normalized score (kcal/mol)", "vina_normalized_score_kcal_mol"),
        ("Vina score basis", "vina_score_basis"), ("Native engine score", "native_score"),
        ("Native score type", "native_score_type"), ("Native details", "native_details_json"),
        ("Native score meanings", "native_score_meanings_json"),
        ("Confidence", "confidence"), ("RMSD to best-energy pose (Å)", "rmsd"),
        ("Atoms compared", "atoms_compared"), ("Pose file in ZIP", "pose_file"),
        ("Run ZIP URL", "bundle_url"), ("Details", "message"),
    ]
    pose_data = []
    for row in poses:
        if not isinstance(row, dict):
            continue
        values = []
        for label, key in pose_headers:
            value = row.get(key)
            if key == "vina_normalized_score_kcal_mol" and value is None:
                value = row.get("common_be", row.get("energy"))
            if key == "vina_score_basis" and value is None:
                value = row.get("common_be_source", row.get("energy_source", row.get("score_type")))
            if key == "native_details_json" and value is None:
                value = json.dumps(row.get("native_details") or {}, sort_keys=True)
            if key == "native_score_meanings_json" and value is None:
                value = score_meanings_json(row.get("engine"), row.get("native_details") or {})
            if key == "rmsd" and isinstance(value, (int, float)):
                value = round(float(value), 3)
            if key == "pose_file" and not value and row.get("pose_index"):
                value = f"poses/pose_{int(row['pose_index']):03d}_normalized.pdb"
            values.append(value)
        pose_data.append(values)

    groups = {}
    for row in poses:
        if not isinstance(row, dict):
            continue
        key = (row.get("batch_id", row.get("run_id")), row.get("target_id"), row.get("ligand_label"), row.get("engine"))
        groups.setdefault(key, []).append(row)
    run_headers = [
        "Target", "Ligand", "Engine", "Run", "Status", "Requested poses", "Returned poses",
        "Best score", "Best pose", "Best-pose RMSD (Å)", "Score type", "Result ZIP URL", "Details",
    ]
    run_data = []
    for key, group in groups.items():
        target, ligand, engine = key[1], key[2], key[3]
        scored = [r for r in group if isinstance(r.get("vina_normalized_score_kcal_mol", r.get("common_be", r.get("energy"))), (int, float))]
        best = min(scored, key=lambda r: float(r.get("vina_normalized_score_kcal_mol", r.get("common_be", r.get("energy"))))) if scored else None
        exemplar = group[0]
        run_data.append([
            target, ligand, exemplar.get("display_name") or engine, key[0], exemplar.get("status"),
            exemplar.get("requested_pose_count"), max((int(r.get("pose_count") or 0) for r in group), default=0),
            (best.get("vina_normalized_score_kcal_mol", best.get("common_be", best.get("energy"))) if best else None),
            best.get("pose_index") if best else None,
            best.get("rmsd", best.get("rmsd_to_best_pose")) if best else None,
            best.get("vina_score_basis", best.get("common_be_source", best.get("energy_source"))) if best else exemplar.get("vina_score_basis", exemplar.get("common_be_source", exemplar.get("energy_source"))),
            next((r.get("bundle_url") for r in group if r.get("bundle_url")), None),
            "; ".join(dict.fromkeys(str(r.get("message")) for r in group if r.get("message"))),
        ])
    pose_widths = [16, 22, 36, 20, 13, 10, 14, 9, 14, 15, 19, 33, 29, 18, 29, 44, 13, 23, 15, 42, 48, 46]
    run_widths = [16, 22, 20, 10, 14, 15, 15, 19, 12, 22, 33, 48, 46]
    return build_workbook([
        {"name": "Pose details", "headers": [x[0] for x in pose_headers], "rows": pose_data, "widths": pose_widths},
        {"name": "Run summary", "headers": run_headers, "rows": run_data, "widths": run_widths},
    ])
