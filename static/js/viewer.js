// ══ CELL 39 · static/js/viewer.js ══
/* ─── viewer.js — 3D viewer (3Dmol): protein, docked pose, grid box ───────────
 *
 *  showStructure(targetId, { pose })   draw protein (+ docked ligand if `pose` given)
 *  drawGridBox(viewer, center, size)   cyan wireframe box with separate X/Y/Z edges
 *  toggleVis('cartoon'|'surface'|'sticks'), toggleSpin(), resetView()
 *
 *  Why the old viewer showed "one green round": the protein was an OPAQUE surface
 *  that hid the ligand sitting in the pocket. Now the protein is a cartoon with a
 *  faint surface, the ligand is drawn as thick sticks, and the pocket residues
 *  around it are shown as sticks too.
 */

const BOX_COLOR = '#00ffcc';
const POCKET_CUTOFF = 5.0;          // Å — residues within this distance of the ligand are shown as sticks
const POSE_INTERACTION_TYPES = {
  hbond:       { color: '#00e5ff', name: 'Hydrogen bond' },
  ionic:       { color: '#ff4d6d', name: 'Salt bridge / charge pair' },
  hydrophobic: { color: '#f4a261', name: 'Hydrophobic contact' },
  pi:          { color: '#b388eb', name: 'π stacking' },
  metal:       { color: '#ffd60a', name: 'Metal coordination' },
};
let viewerLibraryPromise = null;

function clearPoseInteractionLines(render = false) {
  const v = state.viewer;
  (state.interactionShapes || []).forEach(shape => { try { if (v) v.removeShape(shape); } catch (_) {} });
  state.interactionShapes = [];
  if (render && v) v.render();
}

function drawPoseInteractionLines(data, rid) {
  const v = state.viewer;
  if (!v || rid !== state.renderId || !state.interactionsEnabled || !data || !Array.isArray(data.items)) return;
  clearPoseInteractionLines();
  data.items.forEach(item => {
    const meta = POSE_INTERACTION_TYPES[item.type];
    const a = item.lig, b = item.prot;
    if (!meta || !Array.isArray(a) || !Array.isArray(b) || a.length < 3 || b.length < 3) return;
    const start = a.slice(0, 3).map(Number), end = b.slice(0, 3).map(Number);
    if (![...start, ...end].every(Number.isFinite)) return;
    const delta = end.map((value, i) => value - start[i]);
    const distance = Math.hypot(...delta);
    if (!distance) return;
    // Short cylinders separated by gaps are consistently dashed across 3Dmol versions.
    const dash = 0.38, pitch = 0.68;
    for (let offset = 0; offset < distance; offset += pitch) {
      const t0 = offset / distance, t1 = Math.min(1, (offset + dash) / distance);
      const p0 = { x: start[0] + delta[0] * t0, y: start[1] + delta[1] * t0, z: start[2] + delta[2] * t0 };
      const p1 = { x: start[0] + delta[0] * t1, y: start[1] + delta[1] * t1, z: start[2] + delta[2] * t1 };
      state.interactionShapes.push(v.addCylinder({ start: p0, end: p1, radius: 0.045, color: meta.color, fromCap: 1, toCap: 1 }));
    }
  });
  v.render();
}

async function loadPoseInteractions(targetId, pose, rid, viewer, ligandModel) {
  state.poseInteractionsLoading = true;
  state.poseInteractionError = '';
  if (typeof dispLegend === 'function') dispLegend();
  try {
    const data = await apiFetch('/api/insight/interactions', {
      method: 'POST', body: JSON.stringify({ target_id: targetId, docked_pdb: pose }),
    });
    if (rid !== state.renderId || state.viewer !== viewer || state.ligModel !== ligandModel) return;
    state.poseInteractionData = data;
  } catch (error) {
    if (rid !== state.renderId || state.viewer !== viewer) return;
    state.poseInteractionData = null;
    state.poseInteractionError = error.message || 'Could not calculate pose interactions';
  } finally {
    if (rid === state.renderId && state.viewer === viewer) {
      state.poseInteractionsLoading = false;
      if (state.interactionsEnabled && state.poseInteractionData) drawPoseInteractionLines(state.poseInteractionData, rid);
      if (typeof dispLegend === 'function') dispLegend();
    }
  }
}

function togglePoseInteractions(enabled) {
  state.interactionsEnabled = !!enabled;
  if (state.interactionsEnabled && state.poseInteractionData) drawPoseInteractionLines(state.poseInteractionData, state.renderId);
  else clearPoseInteractionLines(true);
  if (typeof dispSync === 'function') dispSync();
  if (typeof dispLegend === 'function') dispLegend();
}

async function loadViewerLibrary() {
  if (window.$3Dmol) return true;
  if (viewerLibraryPromise) return viewerLibraryPromise;
  viewerLibraryPromise = new Promise(resolve => {
    const script = document.createElement('script');
    script.src = 'https://3dmol.org/build/3Dmol-min.js';
    script.async = true;
    script.onload = () => {
      const loaded = !!window.$3Dmol;
      if (!loaded) viewerLibraryPromise = null;
      resolve(loaded);
    };
    script.onerror = () => { viewerLibraryPromise = null; resolve(false); };
    document.head.appendChild(script);
  });
  return viewerLibraryPromise;
}

// ── Grid box (cuboid) ─────────────────────────────────────────────────────────
function drawGridBox(viewer, center, size, dim = false) {
  (state.boxShapes || []).forEach(s => { try { viewer.removeShape(s); } catch (e) {} });
  state.boxShapes = [];
  if (!viewer || !center || !size) return;

  const [cx, cy, cz] = center;
  const [sx, sy, sz] = [size[0], size[1] ?? size[0], size[2] ?? size[0]];
  const hx = sx / 2, hy = sy / 2, hz = sz / 2;

  // (edges only — filled faces hid everything inside the box)

  // 12 solid edges (cylinders are always visible, unlike wireframe flags)
  const corner = (i, j, k) => ({ x: cx + (i ? hx : -hx), y: cy + (j ? hy : -hy), z: cz + (k ? hz : -hz) });
  const edge = (a, b) => state.boxShapes.push(viewer.addCylinder({
    start: a, end: b, radius: dim ? 0.04 : 0.06, color: BOX_COLOR, fromCap: 1, toCap: 1, alpha: 0.6, opacity: 0.6,
  }));
  for (const [j, k] of [[0,0],[0,1],[1,0],[1,1]]) edge(corner(0, j, k), corner(1, j, k));   // along X
  for (const [i, k] of [[0,0],[0,1],[1,0],[1,1]]) edge(corner(i, 0, k), corner(i, 1, k));   // along Y
  for (const [i, j] of [[0,0],[0,1],[1,0],[1,1]]) edge(corner(i, j, 0), corner(i, j, 1));   // along Z
  state.boxModel = true;
}

// ── helpers ───────────────────────────────────────────────────────────────────
function viewerContainer() { return document.getElementById('viewer-container'); }

function clearProteinViewer() {
  state.renderId = (state.renderId || 0) + 1;
  const viewer = state.viewer;
  if (viewer) {
    try {
      viewer.removeAllModels();
      viewer.removeAllShapes();
      viewer.removeAllSurfaces();
      viewer.removeAllLabels();
      viewer.render();
    } catch (_) {}
  }
  state.viewer = null;
  state.proteinModel = null;
  state.ligModel = null;
  state.boxModel = null;
  state.boxShapes = [];
  state.pocketSerials = [];
  state.interactionShapes = [];
  state.hoverSphere = null;
  state.spinning = false;
  const legend = document.getElementById('viewer-legend');
  if (legend) legend.textContent = '';
  const container = viewerContainer();
  if (container) container.innerHTML = `
    <div class="viewer-placeholder">
      <div class="big-icon">🔬</div>
      <div>Select or import a PDB to view its protein structure</div>
      <div style="font-size:11px;color:var(--text-3)">No protein loads until you select a PDB</div>
    </div>
    <div id="coord-overlay" class="coord-overlay" style="display:none">
      X: <span id="coord-x">0.0</span> Y: <span id="coord-y">0.0</span> Z: <span id="coord-z">0.0</span>
    </div>`;
}

function ensureViewer() {
  const container = viewerContainer();
  if (!container) return null;
  if (!window.$3Dmol) {
    container.innerHTML = '<div class="viewer-placeholder"><div class="big-icon">🧬</div><div>3D viewer library could not be loaded.</div></div>';
    return null;
  }
  if (!state.viewer || !container.querySelector('canvas')) {
    container.innerHTML = '';
    state.viewer = $3Dmol.createViewer(container, { backgroundColor: '#020408', antialias: true });
    if (!state.viewer) {
      container.innerHTML = '<div class="viewer-placeholder"><div class="big-icon">🧬</div><div>This browser could not create a WebGL protein viewer.</div></div>';
      return null;
    }

        if (state.viewerRO) { try { state.viewerRO.disconnect(); } catch (e) {} }
    if (window.ResizeObserver) {
      let _t = null;
      state.viewerRO = new ResizeObserver(() => {
        clearTimeout(_t);
        _t = setTimeout(() => {
          if (state.viewer && container.clientWidth > 0) {
            state.viewer.resize();
            state.viewer.render();
          }
        }, 50);
      });
      state.viewerRO.observe(container);
    }

    const overlay = document.createElement('div');
    overlay.id = 'coord-overlay';
    overlay.className = 'coord-overlay';
    overlay.style.display = 'none';
    overlay.innerHTML = 'X: <span id="coord-x">0.0</span> Y: <span id="coord-y">0.0</span> Z: <span id="coord-z">0.0</span>';
    container.appendChild(overlay);

    const legend = document.createElement('div');
    legend.id = 'viewer-legend';
    legend.style.cssText = 'position:absolute;left:10px;bottom:8px;font-size:11px;color:var(--text-2);pointer-events:none;' +
                           'background:rgba(2,4,8,.65);padding:3px 8px;border-radius:6px;z-index:5';
    container.appendChild(legend);
  }
  const canvas = container.querySelector('canvas');
  if (canvas) {
    canvas.style.display = 'block';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
  }
  return state.viewer;
}

async function fetchReceptorPdb(targetId) {
  state.pdbCache = state.pdbCache || {};
  if (state.pdbCache[targetId]) return state.pdbCache[targetId];
  const resp = await fetch(`/api/receptor/${targetId}/pdb`);
  if (!resp.ok) throw new Error(`receptor file not found (${resp.status})`);
  const txt = await resp.text();
  state.pdbCache[targetId] = txt;
  return txt;
}

// Serial numbers of every protein atom in a residue that has an atom within `cutoff` Å of the ligand
function pocketSerials(proteinModel, ligModel, cutoff) {
  const lig = ligModel.selectedAtoms({}).filter(a => a.elem !== 'H');
  const prot = proteinModel.selectedAtoms({});
  const c2 = cutoff * cutoff;
  const hit = new Set();
  for (const p of prot) {
    const key = `${p.chain}:${p.resi}`;
    if (hit.has(key)) continue;
    for (const l of lig) {
      const dx = p.x - l.x, dy = p.y - l.y, dz = p.z - l.z;
      if (dx * dx + dy * dy + dz * dz <= c2) { hit.add(key); break; }
    }
  }
  return prot.filter(p => hit.has(`${p.chain}:${p.resi}`)).map(p => p.serial);
}

function styleLigand(model, thick = true) {
  model.setStyle({}, {
    stick: { colorscheme: 'greenCarbon', radius: thick ? 0.24 : 0.18 },
    sphere: { colorscheme: 'greenCarbon', scale: thick ? 0.28 : 0.2 },
  });
}

function styleProtein(mode = 'both') {
  const v = state.viewer, p = state.proteinModel;
  if (!v || !p) return;
  v.removeAllSurfaces();
  p.setStyle({}, { cartoon: { color: '#7f8fa6', opacity: 0.9, thickness: 0.35 } });
  if (state.pocketSerials && state.pocketSerials.length) {
    p.setStyle({ serial: state.pocketSerials }, {
      cartoon: { color: '#7f8fa6', opacity: 0.9 },
      stick: { colorscheme: 'cyanCarbon', radius: 0.12 },
    });
  }

  if (mode === 'both') {
    try {
      const sp = v.addSurface($3Dmol.SurfaceType.VDW,
        { opacity: 0.16, color: '#aab8c8' }, { model: p.getID() });
      if (sp && sp.catch) sp.catch(err => console.warn('Surface failed:', err));
    } catch (e) { console.warn('Surface failed:', e); }
  }
  // STYLEPROTEIN_FIXED
  v.render();
}

// ── main entry ────────────────────────────────────────────────────────────────
async function showStructure(targetId, { pose = null } = {}) {
  if (!targetId || !state.targets?.[targetId]) return;
  if (!window.$3Dmol) await loadViewerLibrary();
  if (state.selectedTarget !== targetId) return;
  const viewer = ensureViewer();
  if (!viewer) return;
  const container = viewerContainer();

  const rid = state.renderId = (state.renderId || 0) + 1;
  clearPoseInteractionLines();
  state.poseInteractionData = null;
  state.poseInteractionError = '';
  state.poseInteractionsLoading = false;
  if (typeof dispLegend === 'function') dispLegend();
  let receptorPdb;
  try {
    receptorPdb = await fetchReceptorPdb(targetId);
  } catch (e) {
    if (rid === state.renderId) {
      container.innerHTML = `<div class="viewer-placeholder"><div class="big-icon">🧬</div><div>Protein structure could not be loaded.<br><span style="font-size:11px">${escapeHtml(e.message)}</span></div></div>`;
      state.viewer = null;
    }
    showToast(`Could not load the protein for the 3D view: ${e.message}`, 'error', '⚠️');
    return;
  }

  if (rid !== state.renderId) return;

  try {
    viewer.removeAllModels();
    viewer.removeAllShapes();
    viewer.removeAllSurfaces();
    viewer.removeAllLabels();
    state.boxShapes = [];
    state.interactionShapes = [];
    state.hoverSphere = null;
    state.ligModel = null;
    state.pocketSerials = [];

    const protein = viewer.addModel(receptorPdb, 'pdb');
    const proteinAtoms = protein.selectedAtoms({});
    if (!proteinAtoms.length) throw new Error('The receptor file loaded, but it contains no drawable protein atoms.');
    state.proteinModel = protein;

    // docked ligand (always PDB text from the server — all engines are normalised)
    let ligModel = null;
    if (pose && pose.trim()) {
      ligModel = viewer.addModel(pose, 'pdb');
      if (!ligModel.selectedAtoms({}).length) {
        viewer.removeModel(ligModel);
        ligModel = viewer.addModel(pose, 'pdbqt');
      }
      if (!ligModel.selectedAtoms({}).length) {
        viewer.removeModel(ligModel);
        ligModel = null;
        showToast('The docking result did not contain drawable ligand coordinates.', 'error', '⚠️');
      }
    }
    if (ligModel) {
      state.ligModel = ligModel;
      styleLigand(ligModel);
      state.pocketSerials = pocketSerials(protein, ligModel, POCKET_CUTOFF);
    }
    styleProtein('cartoon');

    // grid box
    const t = state.targets[targetId];
    if (t && t.center && t.size) drawGridBox(viewer, t.center, t.size, !!ligModel);

    // click a protein atom → move the grid box centre there
    protein.setClickable({}, true, function (atom) {
      if (!state.selectedTarget) return;
      document.getElementById('box-cx').value = atom.x.toFixed(2);
      document.getElementById('box-cy').value = atom.y.toFixed(2);
      document.getElementById('box-cz').value = atom.z.toFixed(2);
      showToast(`Grid centre → (${atom.x.toFixed(1)}, ${atom.y.toFixed(1)}, ${atom.z.toFixed(1)})`, 'success', '🎯');
      updateBox();
    });

    // hover: coordinate read-out + small marker
    protein.setHoverable({}, true,
      function (atom, v) {
        const o = document.getElementById('coord-overlay');
        if (o) {
          o.style.display = 'block';
          o.querySelector('#coord-x').textContent = atom.x.toFixed(1);
          o.querySelector('#coord-y').textContent = atom.y.toFixed(1);
          o.querySelector('#coord-z').textContent = atom.z.toFixed(1);
        }
        if (state.hoverSphere) v.removeShape(state.hoverSphere);
        state.hoverSphere = v.addSphere({ center: { x: atom.x, y: atom.y, z: atom.z }, radius: 1.0, color: '#4cc9f0', opacity: 0.5 });
        v.render();
      },
      function (atom, v) {
        if (state.hoverSphere) { v.removeShape(state.hoverSphere); state.hoverSphere = null; v.render(); }
        const o = document.getElementById('coord-overlay');
        if (o) o.style.display = 'none';
      });

    // camera
    // Fit all loaded models. A model-only selector can resolve to an empty
    // selection in some 3Dmol builds, leaving a healthy canvas completely black.
    viewer.zoomTo();
    if (ligModel) viewer.zoom(1.15);

    const legend = document.getElementById('viewer-legend');
    if (legend) legend.innerHTML = ligModel
      ? '<span style="color:#2af598">■</span> docked ligand &nbsp; <span style="color:#4cc9f0">■</span> pocket residues &nbsp; <span style="color:#00ffcc">▭</span> grid box'
      : '<span style="color:#00ffcc">▭</span> grid box &nbsp;·&nbsp; click the protein to move its centre';

    viewer.resize();
    viewer.render();
    setTimeout(() => { viewer.resize(); viewer.render(); }, 400);
    if (typeof dispSync === 'function') dispSync();
    if (ligModel && pose && pose.trim()) loadPoseInteractions(targetId, pose, rid, viewer, ligModel);
    else if (typeof dispLegend === 'function') dispLegend();
    if (ligModel) showToast('3D complex loaded — docked ligand shown in green', 'success', '🔬');
  } catch (e) {
    console.error('3Dmol render error:', e);
    container.innerHTML = `<div class="viewer-placeholder"><div class="big-icon">🧬</div>
      <div>3D render error.<br><span style="font-size:11px;color:var(--text-3)">${escapeHtml(e.message)}</span></div></div>`;
    state.viewer = null;
  }
}

// 3Dmol may be hidden while another tab is active. Resize it after the Dock
// panel becomes visible again so its canvas is not left at zero dimensions.
document.addEventListener('click', event => {
  if (!event.target.closest('[data-tab="dock"]')) return;
  requestAnimationFrame(() => {
    if (state.viewer) {
      state.viewer.resize();
      state.viewer.render();
    } else if (state.selectedTarget) {
      showStructure(state.selectedTarget);
    }
  });
});

// ── toolbar buttons ───────────────────────────────────────────────────────────
function previewGridBox(scroll = true) {
  if (!state.selectedTarget) { showToast('Select a protein first', 'error', '⚠️'); return; }
  if (scroll) viewerContainer().parentElement.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  showStructure(state.selectedTarget);
}

function toggleVis(type) {
  if (!state.viewer || !state.proteinModel) return;
  if (type === 'cartoon') styleProtein('cartoon');
  else if (type === 'surface') styleProtein('both');
  else if (type === 'sticks' && state.ligModel) {
    state.sticksThick = !state.sticksThick;
    styleLigand(state.ligModel, state.sticksThick);
  }
  state.viewer.render();
}

function toggleSpin() {
  if (!state.viewer) return;
  state.spinning = !state.spinning;
  state.viewer.spin(state.spinning ? 'y' : false);
  const b = document.getElementById('btn-spin');
  if (b) b.classList.toggle('active', state.spinning);
}

function resetView() {
  if (!state.viewer) return;
  if (state.ligModel) { state.viewer.zoomTo({ model: state.ligModel.getID() }); state.viewer.zoom(0.55); }
  else state.viewer.zoomTo();
  state.viewer.render();
}
