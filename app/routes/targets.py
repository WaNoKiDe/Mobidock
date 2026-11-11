# ══ CELL 28 · app/routes/targets.py ══
"""
routes/targets.py — Targets, grid box, receptor files, protein browser (RCSB), upload, literature.
"""
import asyncio
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import logging
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

from fastapi import APIRouter, File, HTTPException, UploadFile
from fastapi.responses import FileResponse

from .. import docking as dock_engine
from ..config import MAX_BOX_EDGE, RECEPTORS_DIR
from ..receptor_prep import RECEPTORS, RECEPTORS_DIR, _NOT_LIGANDS, download_pdb, find_component_binding_site, process_custom_pdb
from ..runtime import BoxUpdateRequest, PdbImportRequest

logger = logging.getLogger("medchem.targets")
router = APIRouter()
_NATIVE_LIGAND_CACHE = {}
# Native site choices include small organic co-crystal components such as
# tartrate (TLA). Keep excluding common solvents, ions, and crystallization
# additives from receptor preparation, but don't hide them from this selector.
_NATIVE_LIGAND_EXCLUSIONS = _NOT_LIGANDS - {"TLA"}

# ── Targets & receptor files ──────────────────────────────────────────────────

def _target_json(tid, info):
    p = dock_engine.RECEPTOR_PATHS.get(tid)
    box = dock_engine.BOX_PARAMS.get(tid, {})
    return {
        "pdb_id": info.get("pdb_id", tid),
        "name": info.get("name", tid),
        "disease": info.get("disease", ""),
        "color": info.get("color", "#1abc9c"),
        "icon": info.get("icon", "🧪"),
        "center": list(box.get("center", info.get("center", (0, 0, 0)))),
        "size": list(box.get("size", info.get("size", (20, 20, 20)))),
        "site_source": info.get("site_source", "preset"),
        "native_ligand": info.get("native_ligand", ""),
        "site_ligand": info.get("site_ligand", ""),
        "native_ligand_removed": bool(info.get("native_ligand_removed", False)),
        "native_ligands": info.get("native_ligands", []),
        "box_source": info.get("box_source", ""),
        "ready": bool(p and Path(p).exists()),
    }


@router.get("/api/targets")
async def list_targets():
    """All docking targets (built-in + imported) with their current grid box."""
    return {
        tid: _target_json(tid, info)
        for tid, info in RECEPTORS.items()
    }


@router.get("/api/receptor/{target_id}/pdb")
async def get_receptor_pdb(target_id: str):
    """Clean receptor PDB for the 3D viewer."""
    if target_id not in RECEPTORS:
        raise HTTPException(status_code=404, detail=f"Unknown target: {target_id}")
    info = RECEPTORS[target_id]
    pdb_path = Path(info["pdb_path"]) if "pdb_path" in info else RECEPTORS_DIR / f"{target_id}_{info['pdb_id']}_clean.pdb"
    if not pdb_path.exists():
        raise HTTPException(status_code=404, detail=f"Receptor PDB not found for {target_id}.")
    return FileResponse(str(pdb_path), media_type="text/plain", filename=f"{target_id}.pdb")


def _native_ligand_pdb(raw_pdb: str, comp_id: str, chain=None, residue=None, insertion=None) -> str:
    """Extract the chosen co-crystal ligand copy, including explicit CONECT bonds."""
    first_model = (raw_pdb or "").split("ENDMDL", 1)[0]
    groups = {}
    for line in first_model.splitlines():
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        if line[17:20].strip().upper() != str(comp_id or "").strip().upper():
            continue
        key = (line[21].strip() or "_", line[22:26].strip(), line[26].strip())
        if chain is not None and key != (str(chain), str(residue), str(insertion or "")):
            continue
        if len(line) > 16 and line[16] not in (" ", "A", "1"):
            continue
        groups.setdefault(key, []).append(line)
    if not groups:
        return ""
    atoms = max(groups.values(), key=len)
    serials = set()
    for line in atoms:
        try:
            serials.add(int(line[6:11]))
        except ValueError:
            pass
    conect = []
    for line in first_model.splitlines():
        if not line.startswith("CONECT"):
            continue
        try:
            numbers = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
        except ValueError:
            continue
        if numbers and numbers[0] in serials:
            bonded = [number for number in numbers[1:] if number in serials]
            if bonded:
                conect.append("CONECT" + f"{numbers[0]:5d}" + "".join(f"{number:5d}" for number in bonded))
    return "\n".join([*atoms, *conect, "END", ""])


def _ccd_molecule_pdb(mol, center) -> str:
    """Place a CCD reference conformer at the native-site centre for display."""
    try:
        from rdkit.Geometry import Point3D
        if mol is None or mol.GetNumConformers() == 0 or len(center) != 3:
            return ""
        conf = mol.GetConformer()
        heavy = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() != "H"]
        if not heavy:
            return ""
        points = [conf.GetAtomPosition(i) for i in heavy]
        centroid = [sum(point[axis] for point in points) / len(points) for axis in range(3)]
        delta = [float(center[axis]) - centroid[axis] for axis in range(3)]
        for i in range(mol.GetNumAtoms()):
            point = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(*(point[axis] + delta[axis] for axis in range(3))))
        return Chem.MolToPDBBlock(mol)
    except Exception:
        return ""


def _native_ccd_candidate(target_id: str, comp: str) -> dict | None:
    """Build a native-ligand choice from its recorded CCD ID when PDB atoms are absent."""
    comp = (comp or "").strip().upper()
    if not comp:
        return None
    info = RECEPTORS.get(target_id) or {}
    box = dock_engine.BOX_PARAMS.get(target_id) or {}
    center = list(box.get("center", info.get("center", ())))
    size = list(box.get("size", info.get("size", ())))
    if len(center) != 3 or len(size) != 3:
        return None
    try:
        from rdkit import Chem
        url = f"https://files.rcsb.org/ligands/download/{urllib.parse.quote(comp)}.sdf"
        request = urllib.request.Request(url, headers={"User-Agent": "MOBIDock/1.0"})
        with urllib.request.urlopen(request, timeout=8) as response:
            sdf = response.read().decode("utf-8", errors="ignore")
        mol = Chem.MolFromMolBlock(sdf, removeHs=True, sanitize=True)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if not smiles:
            return None
        title = sdf.splitlines()[0].strip() if sdf.splitlines() else ""
        return {
            "ccd_id": comp, "name": title or comp, "smiles": smiles,
            "atom_count": int(Chem.RemoveHs(mol).GetNumAtoms()),
            "chain": "—", "residue": "—", "center": center, "size": size,
            "pdb": _ccd_molecule_pdb(mol, center),
        }
    except Exception as exc:
        logger.debug("Could not load CCD structure %s for %s: %s", comp, target_id, exc)
        return None


def _rcsb_native_component_ids(pdb_id: str) -> list[str]:
    """Get ligand CCD IDs from RCSB when local PDB coordinates are unavailable."""
    pdb_id = (pdb_id or "").strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return []
    query = (
        'query { entry(entry_id: "' + pdb_id + '") { nonpolymer_entities {'
        ' nonpolymer_comp { chem_comp { id formula_weight } } } } }'
    )
    try:
        request = urllib.request.Request(
            "https://data.rcsb.org/graphql",
            data=json.dumps({"query": query}).encode("utf-8"),
            headers={"User-Agent": "MOBIDock/1.0", "Content-Type": "application/json"},
        )
        with urllib.request.urlopen(request, timeout=12) as response:
            payload = json.loads(response.read().decode("utf-8", errors="ignore"))
        entities = ((payload.get("data") or {}).get("entry") or {}).get("nonpolymer_entities") or []
        weighted = []
        for entity in entities:
            chem = ((entity.get("nonpolymer_comp") or {}).get("chem_comp") or {})
            comp = str(chem.get("id") or "").strip().upper()
            if not comp or comp in _NATIVE_LIGAND_EXCLUSIONS:
                continue
            try:
                weight = float(chem.get("formula_weight") or 0)
            except (TypeError, ValueError):
                weight = 0.0
            if weight >= 80:
                weighted.append((weight, comp))
        if weighted:
            return [comp for _, comp in sorted(set(weighted), reverse=True)]
    except Exception as exc:
        logger.info("Could not read native-ligand GraphQL metadata for PDB %s: %s", pdb_id, exc)

    # Keep a REST fallback: the app already uses these stable core endpoints for
    # the PDB browser, and some RCSB deployments reject the GraphQL query above.
    try:
        entry = _related_rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=12)
        entity_ids = (
            entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])
            or []
        )
        components = []
        for entity_id in entity_ids:
            try:
                entity = _related_rcsb_json(
                    f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}",
                    timeout=8,
                )
                component = (entity.get("pdbx_entity_nonpoly") or {}).get("comp_id")
                comp = str(component or "").strip().upper()
                if comp and comp not in _NATIVE_LIGAND_EXCLUSIONS and comp not in components:
                    components.append(comp)
            except Exception as exc:
                logger.debug("Could not read ligand metadata for %s entity %s: %s", pdb_id, entity_id, exc)
        return components
    except Exception as exc:
        logger.info("Could not read native-ligand metadata for PDB %s: %s", pdb_id, exc)
        return []


def _native_ligand_candidates(target_id: str) -> list[dict]:
    """Find dockable co-crystal compounds and their CCD/PDB-derived structures."""
    if target_id not in RECEPTORS:
        raise ValueError(f"Unknown target: {target_id}")
    info = RECEPTORS[target_id]
    cached = _NATIVE_LIGAND_CACHE.get(target_id)
    if cached is not None:
        return cached

    raw_path = Path(info.get("raw_pdb_path") or RECEPTORS_DIR / f"{target_id}_source.pdb")
    raw = raw_path.read_text(encoding="utf-8", errors="ignore") if raw_path.exists() else ""
    saved_native = None
    if not raw:
        try:
            from .. import recorddb
            saved_native = recorddb.get_native(target_id)
            raw = str((saved_native or {}).get("pdb") or "")
        except Exception as exc:
            logger.debug("Could not read saved native ligand for %s: %s", target_id, exc)
    pdb_id = str(info.get("pdb_id", ""))
    if not raw and len(pdb_id) == 4 and pdb_id.isalnum():
        raw = download_pdb(pdb_id) or ""
    if not raw:
        preferred = str(info.get("native_ligand") or info.get("site_ligand") or
                        (saved_native or {}).get("comp") or "").upper()
        component_ids = ([preferred] if preferred else []) + _rcsb_native_component_ids(pdb_id)
        candidates = []
        for comp in dict.fromkeys(component_ids):
            fallback = _native_ccd_candidate(target_id, comp)
            if fallback:
                candidates.append(fallback)
        if candidates:
            _NATIVE_LIGAND_CACHE[target_id] = candidates
            return candidates
        return []
    # Keep the source PDB next to the prepared receptor. Built-in receptors are
    # normally stored ligand-free, so this cached source is what makes native
    # co-crystal ligands available offline after the first successful lookup.
    if not raw_path.exists():
        try:
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            raw_path.write_text(raw, encoding="utf-8")
        except OSError as exc:
            logger.debug("Could not cache source PDB for %s: %s", target_id, exc)

    # PDB files with multiple NMR models should expose the first model only,
    # matching the receptor coordinates that are actually docked against.
    first_model = raw.split("ENDMDL", 1)[0]
    groups = {}
    for line in first_model.splitlines():
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        if len(line) > 16 and line[16] not in (" ", "A", "1"):
            continue
        comp = line[17:20].strip().upper()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = "".join(ch for ch in line[12:16].upper() if ch.isalpha())[:1]
        if not comp or comp in _NATIVE_LIGAND_EXCLUSIONS or element in {"H", "D"}:
            continue
        try:
            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        except ValueError:
            continue
        key = (comp, line[21].strip() or "_", line[22:26].strip(), line[26].strip())
        row = groups.setdefault(key, {"atoms": [], "carbon": 0})
        row["atoms"].append(xyz)
        row["carbon"] += element == "C"

    from rdkit import Chem
    import urllib.error
    import urllib.request

    candidates = []
    # Keep each co-crystal residue as a separate choice. A component can occur
    # in multiple chains or binding sites, and each copy needs its own grid box.
    instances = []
    for key, row in groups.items():
        comp, chain, residue, insertion = key
        if len(row["atoms"]) < 5 or not row["carbon"]:
            continue
        instances.append((comp, {**row, "chain": chain, "residue": residue, "insertion": insertion}))

    preferred = str(info.get("native_ligand") or info.get("site_ligand") or
                    (saved_native or {}).get("comp") or "").upper()
    instances.sort(key=lambda entry: (
        entry[0] != preferred, -len(entry[1]["atoms"]), entry[0],
        entry[1]["chain"], entry[1]["residue"], entry[1]["insertion"],
    ))
    ccd_structures = {}
    for comp, item in instances:
        if comp not in ccd_structures:
            mol, name = None, comp
            try:
                url = f"https://files.rcsb.org/ligands/download/{urllib.parse.quote(comp)}.sdf"
                request = urllib.request.Request(url, headers={"User-Agent": "MOBIDock/1.0"})
                with urllib.request.urlopen(request, timeout=5) as response:
                    sdf = response.read().decode("utf-8", errors="ignore")
                mol = Chem.MolFromMolBlock(sdf, removeHs=True, sanitize=True)
                title = sdf.splitlines()[0].strip() if sdf.splitlines() else ""
                if title:
                    name = title
            except Exception:
                mol = None
            ccd_structures[comp] = (mol, name)
        mol, name = ccd_structures[comp]

        # PDB uploads may contain an unregistered component; infer bonds from the
        # isolated residue coordinates as a local fallback.
        if mol is None:
            residue_lines, serials = [], set()
            for line in first_model.splitlines():
                if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
                    key = (line[17:20].strip().upper(), line[21].strip() or "_", line[22:26].strip(), line[26].strip())
                    if key == (comp, item["chain"], item["residue"], item["insertion"]):
                        alt = line[16] if len(line) > 16 else " "
                        if alt in (" ", "A", "1"):
                            residue_lines.append(line)
                            try: serials.add(int(line[6:11]))
                            except ValueError: pass
            for line in first_model.splitlines():
                if not line.startswith("CONECT"):
                    continue
                try:
                    values = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
                except ValueError:
                    continue
                if len(values) > 1 and values[0] in serials:
                    connected = [v for v in values[1:] if v in serials]
                    if connected:
                        residue_lines.append("CONECT" + f"{values[0]:5d}" + "".join(f"{v:5d}" for v in connected))
            if residue_lines:
                try:
                    mol = Chem.MolFromPDBBlock("\n".join(residue_lines) + "\nEND\n", sanitize=False, removeHs=True, proximityBonding=True)
                    if mol is not None:
                        Chem.SanitizeMol(mol)
                except Exception:
                    mol = None
        if mol is None:
            continue
        try:
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True)
            if not smiles:
                continue
            site = find_component_binding_site(
                first_model, comp, chain=item["chain"], residue=item["residue"], insertion=item["insertion"])
            native_pdb = _native_ligand_pdb(first_model, comp, item["chain"], item["residue"], item["insertion"])
            if not native_pdb:
                continue
            candidates.append({
                "ccd_id": comp, "name": name, "smiles": smiles,
                "atom_count": int(Chem.RemoveHs(mol).GetNumAtoms()),
                "chain": item["chain"], "residue": item["residue"], "insertion": item["insertion"],
                "center": list(site["center"]), "size": list(site["size"]),
                "pdb": native_pdb,
            })
        except Exception as exc:
            logger.debug("Could not build native ligand %s for %s: %s", comp, target_id, exc)

    if not candidates or (preferred and not any(item["ccd_id"] == preferred for item in candidates)):
        fallback_ids = ([preferred] if preferred else []) + _rcsb_native_component_ids(pdb_id)
        for comp in dict.fromkeys(fallback_ids):
            if any(item["ccd_id"] == comp for item in candidates):
                continue
            fallback = _native_ccd_candidate(target_id, comp)
            if fallback:
                candidates.append(fallback)
    if candidates:
        _NATIVE_LIGAND_CACHE[target_id] = candidates
    return candidates


@router.get("/api/targets/{target_id}/native-ligands")
async def native_ligands(target_id: str):
    try:
        return {"target_id": target_id, "ligands": await asyncio.to_thread(_native_ligand_candidates, target_id)}
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@router.post("/api/target/{target_id}/box")
async def update_target_box(target_id: str, req: BoxUpdateRequest):
    """Set the grid box (centre + X/Y/Z edge lengths in Å) for a target."""
    if target_id not in RECEPTORS:
        raise HTTPException(status_code=404, detail=f"Target {target_id} not found")
    if len(req.center) != 3 or len(req.size) != 3:
        raise HTTPException(status_code=422, detail="center and size must each have 3 values")
    if any(v < 1.0 for v in req.size) or any(v > MAX_BOX_EDGE for v in req.size):
        raise HTTPException(status_code=422, detail=f"Grid box edges must be between 1 and {MAX_BOX_EDGE:.0f} Å")
    c, s = dock_engine.set_box(target_id, req.center, req.size)
    RECEPTORS[target_id]["center"], RECEPTORS[target_id]["size"] = c, s
    return {"status": "ok", "center": list(c), "size": list(s)}


# ── Protein browser (RCSB PDB) ────────────────────────────────────────────────

@router.get("/api/pdb/search")
async def search_pdb(q: str, page: int = 1, page_size: int = 10):
    """Search RCSB PDB by keyword (RCSB Search API v2, full_text service) + details via GraphQL."""
    q = (q or "").strip()
    if not q:
        raise HTTPException(status_code=400, detail="Query is required")
    page, page_size = max(1, page), max(1, min(page_size, 25))
    start = (page - 1) * page_size

    def post_json(url, payload, timeout=20):
        rq = urllib.request.Request(
            url, data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "MedChemGame/3.1"})
        with urllib.request.urlopen(rq, timeout=timeout) as resp:
            body = resp.read().decode("utf-8", errors="ignore").strip()
            return json.loads(body) if body else {}

    def run_search():
        node = {"type": "terminal", "service": "full_text", "parameters": {"value": q}}
        variants = [
            {"paginate": {"start": start, "rows": page_size}, "results_content_type": ["experimental"],
             "sort": [{"sort_by": "score", "direction": "desc"}]},
            {"paginate": {"start": start, "rows": page_size}},          # plain fallback
        ]
        last = None
        for opts in variants:
            try:
                return post_json("https://search.rcsb.org/rcsbsearch/v2/query",
                                 {"query": node, "return_type": "entry", "request_options": opts})
            except urllib.error.HTTPError as e:
                if e.code == 204:
                    return {}
                last = e
                try: detail = e.read().decode("utf-8", errors="ignore")[:200]
                except Exception: detail = ""
                logger.warning("RCSB search HTTP %s: %s", e.code, detail)
        raise last

    def details(ids):
        gql = ("query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} "
               "rcsb_entry_info{resolution_combined experimental_method} "
               "polymer_entities{rcsb_entity_source_organism{scientific_name}}}}")
        data = post_json("https://data.rcsb.org/graphql", {"query": gql, "variables": {"ids": ids}}, timeout=15)
        by_id = {}
        for e in ((data.get("data") or {}).get("entries") or []):
            if not e: continue
            info = e.get("rcsb_entry_info") or {}
            res = info.get("resolution_combined") or [None]
            organism = "Unknown"
            for ent in (e.get("polymer_entities") or []):
                src = (ent or {}).get("rcsb_entity_source_organism") or []
                if src and src[0].get("scientific_name"):
                    organism = src[0]["scientific_name"]; break
            by_id[e.get("rcsb_id", "")] = {
                "pdb_id": e.get("rcsb_id", ""),
                "title": (e.get("struct") or {}).get("title") or "Unknown",
                "organism": organism,
                "resolution": round(res[0], 2) if res[0] else None,
                "method": info.get("experimental_method") or "",
            }
        return by_id

    def work():
        data = run_search()
        ids = [h["identifier"] for h in data.get("result_set", [])]
        total = data.get("total_count", 0)
        if not ids:
            return {"query": q, "total": 0, "page": page, "results": []}
        try:
            by_id = details(ids)
        except Exception as exc:
            logger.warning("RCSB GraphQL detail fetch failed: %s", exc)
            by_id = {}
        results = [by_id.get(i) or {"pdb_id": i, "title": "", "organism": "", "resolution": None, "method": ""} for i in ids]
        return {"query": q, "total": total, "page": page, "results": results}

    try:
        return await asyncio.to_thread(work)
    except urllib.error.HTTPError as e:
        logger.exception("RCSB search failed")
        raise HTTPException(status_code=502, detail=f"RCSB API error: {e.code}")
    except Exception as e:
        logger.exception("PDB search failed")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/pdb/info/{pdb_id}")
async def get_pdb_info(pdb_id: str):
    """Return RCSB metadata used in the selected-structure details panel."""
    pdb_id = pdb_id.strip().upper()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        req_obj = urllib.request.Request(url, headers={"User-Agent": "MOBIDock/1.0"})
        with urllib.request.urlopen(req_obj, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="ignore"))
        info = data.get("rcsb_entry_info", {})
        organism = ""
        try:
            query = ("query($ids:[String!]!){entries(entry_ids:$ids){polymer_entities{"
                     "rcsb_entity_source_organism{scientific_name}}}}")
            gql_req = urllib.request.Request(
                "https://data.rcsb.org/graphql",
                data=json.dumps({"query": query, "variables": {"ids": [pdb_id]}}).encode("utf-8"),
                headers={"Content-Type": "application/json", "Accept": "application/json",
                         "User-Agent": "MOBIDock/1.0"},
            )
            with urllib.request.urlopen(gql_req, timeout=8) as resp:
                gql_data = json.loads(resp.read().decode("utf-8", errors="ignore"))
            entries = ((gql_data.get("data") or {}).get("entries") or [])
            for entity in (entries[0] or {}).get("polymer_entities", []) if entries else []:
                sources = (entity or {}).get("rcsb_entity_source_organism") or []
                if sources and sources[0].get("scientific_name"):
                    organism = sources[0]["scientific_name"]
                    break
        except Exception:
            logger.debug("Could not fetch source organism for %s", pdb_id)
        binding_ligands = []
        try:
            summary = await asyncio.to_thread(_pdb_summary_and_ligands, pdb_id)
            binding_ligands = summary.get("binding_ligands", [])
        except Exception:
            logger.debug("Could not fetch bound-ligand details for %s", pdb_id)
        return {
            "pdb_id": pdb_id,
            "title": data.get("struct", {}).get("title", ""),
            "method": data.get("exptl", [{}])[0].get("method", ""),
            "resolution": info.get("resolution_combined", [None])[0],
            "deposition_date": data.get("rcsb_accession_info", {}).get("deposit_date", ""),
            "keywords": data.get("struct_keywords", {}).get("pdbx_keywords", ""),
            "organism": organism,
            "polymer_entity_count": info.get("polymer_entity_count"),
            "non_polymer_entity_count": info.get("nonpolymer_entity_count"),
            "non_polymer_entities": data.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", []),
            "binding_ligands": binding_ligands,
        }
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise HTTPException(status_code=404, detail=f"PDB ID {pdb_id} not found")
        raise HTTPException(status_code=502, detail=f"RCSB API error: {exc.code}")
    except Exception as exc:
        raise HTTPException(status_code=502, detail=str(exc))


def _register(target_id: str, raw_pdb: str, name: str, disease: str) -> dict:
    """Prepare a PDB, add it to RECEPTORS and dock_engine, return its info."""
    info = process_custom_pdb(target_id, raw_pdb)
    info["name"], info["disease"] = name[:80], disease
    RECEPTORS[target_id] = info
    dock_engine.register_custom_target(target_id, info["pdbqt_path"], info["center"], info["size"])
    return info


@router.post("/api/pdb/import")
async def import_pdb(req: PdbImportRequest):
    """Import a PDB entry from RCSB, prepare it for docking and auto-place the grid box."""
    pdb_id = req.pdb_id.strip()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise HTTPException(status_code=400, detail="Invalid PDB ID format (expected 4 alphanumeric characters)")

    def work():
        raw = download_pdb(pdb_id)
        if not raw:
            raise HTTPException(status_code=404, detail=f"PDB ID {pdb_id.upper()} could not be downloaded from RCSB")
        title = f"Protein {pdb_id.upper()}"
        try:
            rq = urllib.request.Request(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}", headers={"User-Agent": "MedChemGame/3.1"})
            with urllib.request.urlopen(rq, timeout=6) as r:
                title = json.loads(r.read().decode()).get("struct", {}).get("title", title)
        except Exception:
            pass
        tid = f"PDB_{pdb_id.upper()}_{uuid.uuid4().hex[:4].upper()}"
        return tid, _register(tid, raw, title, "Imported from RCSB PDB")

    try:
        tid, info = await asyncio.to_thread(work)
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("Failed to process downloaded PDB")
        raise HTTPException(status_code=500, detail=str(e))
    info["pdb_id"] = pdb_id.upper()
    return {"status": "ok", "target_id": tid, "info": _target_json(tid, info)}


@router.post("/api/upload-receptor")
async def upload_receptor(file: UploadFile = File(...)):
    """Upload a custom PDB file as a new target."""
    if not (file.filename or "").lower().endswith(".pdb"):
        raise HTTPException(status_code=400, detail="Only .pdb files are supported")
    raw = (await file.read()).decode("utf-8", errors="ignore")
    tid = f"CUSTOM_{uuid.uuid4().hex[:6].upper()}"
    try:
        info = await asyncio.to_thread(_register, tid, raw, "Custom Target", "User Upload")
    except Exception as e:
        logger.exception("Failed to process uploaded PDB")
        raise HTTPException(status_code=500, detail=str(e))
    return {"status": "ok", "target_id": tid, "info": _target_json(tid, info)}


# ── Literature ────────────────────────────────────────────────────────────────

@router.get("/api/literature")
async def get_literature(target_id: str):
    """Fetch literature related to the target from Europe PMC."""
    if target_id not in RECEPTORS:
        raise HTTPException(status_code=404, detail="Target not found")

    target_name = RECEPTORS[target_id].get("name", target_id)
    query = f'("{target_name}" OR "{target_id}") AND (inhibitor OR "binding affinity" OR IC50)'
    encoded_query = urllib.parse.quote(query)

    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={encoded_query}&format=json&resultType=lite&pageSize=5"

    try:
        req_obj = urllib.request.Request(url, headers={"User-Agent": "MedChemGame/2.0"})
        with urllib.request.urlopen(req_obj) as response:
            data = json.loads(response.read().decode())
            results = data.get("resultList", {}).get("result", [])
            papers = []
            for r in results:
                papers.append({
                    "title": r.get("title"),
                    "authors": r.get("authorString", "Unknown"),
                    "journal": r.get("journalTitle", "Unknown Journal"),
                    "year": r.get("pubYear", ""),
                    "doi": r.get("doi", ""),
                    "pmid": r.get("pmid", ""),
                })
            return {"target": target_name, "papers": papers}
    except Exception as e:
        logger.exception("Europe PMC query failed")
        return {"target": target_name, "papers": [], "error": str(e)}


# ── Routes: Leaderboard ──────────────────────────────────────────────────────


# === LIGAND_PDB_RCSB_V1 ===
# Exact ligand-to-PDB search. Uses RCSB's chemical Search API and Data API.
import time
from typing import Any, Dict, List
from pydantic import BaseModel
from rdkit import Chem

_LIGAND_PDB_SEARCH_CACHE: Dict[str, tuple] = {}
_LIGAND_PDB_TARGET_CACHE: Dict[tuple, str] = {}


class LigandPdbImportRequest(BaseModel):
    pdb_id: str
    ccd_id: str
    smiles: str



def _rcsb_json(url: str, timeout: int = 20) -> dict:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "MedChemGame/3.1 ligand-PDB browser"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        return _rcsb_read_json_safe(response)


def _canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise ValueError("The loaded SMILES is not valid.")
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _search_rcsb_for_exact_ligand(smiles: str, limit: int = 6) -> dict:
    """
    Search only for PDB entries containing the exact ligand structure.
    No PubChem/ChEMBL results are used here.
    """
    canonical = _canonical_smiles(smiles)
    limit = max(1, min(6, int(limit)))

    cached = _LIGAND_PDB_SEARCH_CACHE.get(canonical)
    if cached and time.time() - cached[0] < 600:
        return cached[1]

    query = {
        "query": {
            "type": "terminal",
            "service": "chemical",
            "parameters": {
                "value": canonical,
                "type": "descriptor",
                "descriptor_type": "SMILES",
                "match_type": "graph-exact",
            },
        },
        "return_type": "non_polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 20},
            "results_content_type": ["experimental"],
            "scoring_strategy": "chemical",
        },
    }

    request = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(query).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "MedChemGame/3.1 ligand-PDB browser",
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=12) as response:
            search_data = _rcsb_read_json_safe(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 204:
            answer = {"query_smiles": canonical, "results": []}
            _LIGAND_PDB_SEARCH_CACHE[canonical] = (time.time(), answer)
            return answer
        raise

    candidates = []
    seen_pdb_ids = set()
    # Limit detail lookups to unique structures and fetch them concurrently.
    # RCSB detail calls were previously serial, leaving the UI busy for minutes
    # when one of the returned structures was slow or unavailable.
    for hit in search_data.get("result_set", []):
        identifier = str(hit.get("identifier", ""))
        if "_" not in identifier:
            continue
        pdb_id, entity_id = identifier.rsplit("_", 1)
        pdb_id = pdb_id.upper()
        if pdb_id in seen_pdb_ids or len(pdb_id) != 4:
            continue
        seen_pdb_ids.add(pdb_id)
        candidates.append((pdb_id, entity_id))
        if len(candidates) >= max(12, limit * 2):
            break

    def describe(hit):
        pdb_id, entity_id = hit
        ligand_info = _rcsb_json(
            f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}",
            timeout=4,
        )
        component = ligand_info.get("pdbx_entity_nonpoly", {})
        ccd_id = str(component.get("comp_id", "")).strip().upper()
        if not ccd_id:
            return None
        entry = _rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=4)
        entry_info = entry.get("rcsb_entry_info", {})
        resolutions = entry_info.get("resolution_combined", []) or []
        resolution = resolutions[0] if resolutions else None
        exptl = entry.get("exptl", []) or []
        return {
            "pdb_id": pdb_id,
            "entity_id": str(entity_id),
            "ccd_id": ccd_id,
            "ligand_name": component.get("name") or ccd_id,
            "title": entry.get("struct", {}).get("title") or f"PDB {pdb_id}",
            "resolution": round(float(resolution), 2) if resolution else None,
            "method": exptl[0].get("method", "") if exptl else "",
            "match": "exact co-crystal ligand",
        }

    results = []
    with ThreadPoolExecutor(max_workers=max(1, min(12, len(candidates)))) as pool:
        futures = {pool.submit(describe, hit): hit for hit in candidates}
        for future in as_completed(futures):
            try:
                row = future.result()
                if row:
                    results.append(row)
            except Exception as exc:
                logger.warning("Skipping RCSB result %s_%s: %s", *futures[future], exc)

    # Lower resolution value is preferable.
    results.sort(key=lambda row: (
        row["resolution"] is None,
        row["resolution"] if row["resolution"] is not None else 999.0,
        row["pdb_id"],
    ))

    answer = {"query_smiles": canonical, "results": results[:limit]}
    _LIGAND_PDB_SEARCH_CACHE[canonical] = (time.time(), answer)
    return answer


def _prepare_ligand_selected_target(pdb_id: str, ccd_id: str, smiles: str):
    pdb_id = (pdb_id or "").strip().upper()
    ccd_id = (ccd_id or "").strip().upper()
    canonical = _canonical_smiles(smiles)

    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError("Invalid PDB ID.")
    if not ccd_id:
        raise ValueError("The selected PDB result has no co-crystal ligand ID.")

    cache_key = (canonical, pdb_id, ccd_id)
    cached_target = _LIGAND_PDB_TARGET_CACHE.get(cache_key)
    if cached_target and cached_target in RECEPTORS:
        return cached_target, RECEPTORS[cached_target]

    raw_pdb = download_pdb(pdb_id)
    if not raw_pdb:
        raise ValueError(f"RCSB could not download PDB {pdb_id}.")

    target_id = f"LIGAND_{pdb_id}_{uuid.uuid4().hex[:6].upper()}"

    # process_custom_pdb cleans all HETATM records: the native ligand is removed.
    # ligand_comp_id separately preserves its coordinates only for grid placement.
    info = process_custom_pdb(target_id, raw_pdb, ligand_comp_id=ccd_id)
    info["pdb_id"] = pdb_id
    info["name"] = f"{pdb_id} · native ligand {ccd_id}"
    info["disease"] = "RCSB ligand-selected structure"
    info["native_ligand"] = ccd_id
    info["native_ligand_removed"] = True
    ligand_mol = Chem.MolFromSmiles(canonical)
    info["native_ligands"] = [{
        "ccd_id": ccd_id,
        "name": ccd_id,
        "smiles": canonical,
        "atom_count": ligand_mol.GetNumHeavyAtoms(),
        "chain": "—",
        "residue": "—",
        "center": list(info["center"]),
        "size": list(info["size"]),
        "pdb": _native_ligand_pdb(raw_pdb, ccd_id),
    }]

    RECEPTORS[target_id] = info
    dock_engine.register_custom_target(
        target_id,
        info["pdbqt_path"],
        info["center"],
        info["size"],
    )
    _LIGAND_PDB_TARGET_CACHE[cache_key] = target_id
    return target_id, info


@router.get("/api/ligand-pdb/search")
async def ligand_pdb_search(smiles: str, limit: int = 6):
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(_search_rcsb_for_exact_ligand, smiles, limit),
            timeout=22,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except asyncio.TimeoutError as exc:
        raise HTTPException(status_code=504, detail="RCSB ligand search timed out. Retry the search in a moment.") from exc
    except Exception as exc:
        logger.exception("RCSB ligand search failed")
        raise HTTPException(status_code=502, detail=f"RCSB ligand search failed: {exc}")


@router.post("/api/ligand-pdb/import")
async def ligand_pdb_import(req: LigandPdbImportRequest):
    try:
        target_id, info = await asyncio.to_thread(
            _prepare_ligand_selected_target,
            req.pdb_id,
            req.ccd_id,
            req.smiles,
        )
        return {
            "status": "ok",
            "target_id": target_id,
            "info": _target_json(target_id, info),
            "native_ligand": req.ccd_id.upper(),
            "native_ligand_removed": True,
        }
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception as exc:
        logger.exception("Ligand-selected PDB import failed")
        raise HTTPException(status_code=500, detail=str(exc))


# === END LIGAND_PDB_RCSB_V1 ===


# BEGIN RELATED_PDBS_V1
def _related_rcsb_json(url, timeout=20):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "MedChemGame/3.1", "Accept": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return _rcsb_read_json_safe(response)


def _pdb_summary_and_ligands(pdb_id):
    entry = _related_rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    entry_info = entry.get("rcsb_entry_info", {})
    resolutions = entry_info.get("resolution_combined", []) or []
    exptl = entry.get("exptl", []) or []

    ligand_ids = (
        entry.get("rcsb_entry_container_identifiers", {})
        .get("non_polymer_entity_ids", [])
        or []
    )

    ligands = []
    for entity_id in ligand_ids[:8]:
        try:
            entity = _related_rcsb_json(
                f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}",
                timeout=12,
            )
            component = entity.get("pdbx_entity_nonpoly", {})
            ccd_id = component.get("comp_id")
            if ccd_id:
                ligands.append({
                    "ccd_id": ccd_id,
                    "name": component.get("name") or ccd_id,
                })
        except Exception:
            continue

    return {
        "pdb_id": pdb_id,
        "title": entry.get("struct", {}).get("title") or f"PDB {pdb_id}",
        "method": exptl[0].get("method", "") if exptl else "",
        "resolution": round(float(resolutions[0]), 2) if resolutions else None,
        "binding_ligands": ligands,
    }


def _related_pdbs(pdb_id, limit=6):
    pdb_id = (pdb_id or "").strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError("A PDB ID must contain four letters or numbers.")

    # RCSB structural-similarity search using assembly 1 of the selected PDB.
    payload = {
        "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
                "value": {"entry_id": pdb_id, "assembly_id": "1"},
                "number_of_candidates": 250,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max(12, int(limit) + 2)},
            "results_content_type": ["experimental"],
            "scoring_strategy": "structure",
        },
    }

    request = urllib.request.Request(
        "https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": "MedChemGame/3.1",
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=35) as response:
            data = _rcsb_read_json_safe(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 204:
            return {"pdb_id": pdb_id, "results": []}
        raise

    results = []
    for hit in data.get("result_set", []):
        related_id = str(hit.get("identifier", "")).split("-")[0].upper()
        if len(related_id) != 4 or related_id == pdb_id:
            continue

        try:
            row = _pdb_summary_and_ligands(related_id)
            row["structure_score"] = hit.get("score")
            results.append(row)
        except Exception:
            continue

        if len(results) >= max(1, min(int(limit), 8)):
            break

    return {"pdb_id": pdb_id, "results": results}


@router.get("/api/pdb/{pdb_id}/related")
async def pdb_related_structures(pdb_id: str, limit: int = 6):
    try:
        return await asyncio.to_thread(_related_pdbs, pdb_id, limit)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except Exception as exc:
        logger.warning("Related PDB search failed for %s: %s", pdb_id, exc)
        return {
            "pdb_id": pdb_id.upper(),
            "results": [],
            "warning": "RCSB structural-similarity search is temporarily unavailable.",
        }
# END RELATED_PDBS_V1


# RCSB_SAFE_JSON_V1
def _rcsb_read_json_safe(response):
    import json

    # A successful search with no matching structures has no JSON body.
    if response.status == 204:
        return {"result_set": []}

    raw = response.read().decode("utf-8-sig").strip()
    if not raw:
        raise ValueError(
            "RCSB returned an unexpected empty response. Please retry."
        )

    try:
        result = json.loads(raw)
    except json.JSONDecodeError:
        raise ValueError(
            "RCSB returned a non-JSON response. Please retry."
        ) from None

    if not isinstance(result, dict):
        raise ValueError("RCSB returned an unexpected response format.")

    return result
