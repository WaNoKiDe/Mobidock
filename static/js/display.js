/* ─── display.js — protein display options, colour legend, hydrogens, safe surfaces, ligand preview ─── */

const DISP = {
  repr: 'cartoon', colorBy: 'residue', sideChains: 'auto', surface: 'off', opacity: 0.5,
  hydrogens: false, labels: false, showProtein: true, showLigand: true, showBox: true, preview: true,
  labelObjs: [], surfaceIds: [], tok: 0, previewData: null, previewPending: '', previewModel: null, zoomedFor: '', timer: null,
};

const DISP_CLASSES = [
  { key: 'acidic', name: 'Acidic (−)', res: ['ASP', 'GLU'], color: '#ff4d6d', note: 'Asp, Glu — negative charge' },
  { key: 'basic', name: 'Basic (+)', res: ['LYS', 'ARG'], color: '#4895ef', note: 'Lys, Arg — positive charge' },
  { key: 'polar', name: 'Polar', res: ['SER', 'THR', 'ASN', 'GLN', 'CYS', 'HIS'], color: '#2ec4b6', note: 'Ser, Thr, Asn, Gln, Cys, His — H-bond partners' },
  { key: 'hydrophobic', name: 'Greasy', res: ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO'], color: '#f4a261', note: 'Ala, Val, Leu, Ile, Met, Pro — oily' },
  { key: 'aromatic', name: 'Aromatic', res: ['PHE', 'TRP', 'TYR'], color: '#b388eb', note: 'Phe, Trp, Tyr — flat rings' },
  { key: 'small', name: 'Glycine', res: ['GLY'], color: '#8d99ae', note: 'tiny and flexible' },
];
const DISP_SURFACE_RESIDUE_COLORS = Object.fromEntries(
  DISP_CLASSES.flatMap(group => group.res.map(resn => [resn, group.color]))
);
const DISP_ELEMENTS = [['Carbon', '#909090'], ['Nitrogen', '#3050F8'], ['Oxygen', '#FF0D0D'], ['Sulfur', '#FFFF30'], ['Hydrogen', '#FFFFFF'], ['Iron / metals', '#E06633']];
const DISP_BACKBONE = ['N', 'CA', 'C', 'O', 'OXT'];

const dh = s => escapeHtml(String(s === undefined || s === null ? '' : s));
const dispHexInt = c => parseInt(String(c).replace('#', ''), 16);
function dispResColor(resn) { for (const c of DISP_CLASSES) if (c.res.includes(resn)) return c.color; return null; }

// ── colour/representation helpers ────────────────────────────────────────────
function dispStyleObj(spec) {
  if (DISP.repr === 'cartoon') return { cartoon: Object.assign({ opacity: 0.9, thickness: 0.35 }, spec) };
  if (DISP.repr === 'sticks') return { stick: Object.assign({ radius: 0.14 }, spec) };
  if (DISP.repr === 'lines') return { line: Object.assign({}, spec) };
  return { sphere: Object.assign({ scale: 0.55 }, spec) };
}
function dispBaseSpec() {
  const cart = DISP.repr === 'cartoon';
  if (DISP.colorBy === 'element') return cart ? { color: '#9aa5b1' } : { colorscheme: 'Jmol' };
  if (DISP.colorBy === 'chain') return { colorscheme: 'chain' };
  if (DISP.colorBy === 'spectrum') return { color: 'spectrum' };
  return { color: '#9fb0c3' };
}
function dispStickSpec() {
  if (DISP.colorBy === 'element') return { colorscheme: 'Jmol' };
  if (DISP.colorBy === 'chain') return { colorscheme: 'chain' };
  if (DISP.colorBy === 'spectrum') return { color: 'spectrum' };
  return { color: '#9fb0c3' };
}

// residues lining the grid box (computed in the browser from the loaded protein — no server call needed)
function dispPocketGroups() {
  const t = state.targets && state.targets[state.selectedTarget], p = state.proteinModel;
  if (!t || !t.center || !t.size || !p) return [];
  const c = t.center, h = t.size.map(s => s / 2 + 2.0), seen = {};
  p.selectedAtoms({}).forEach(a => {
    if (a.elem === 'H' || a.hetflag) return;
    if (Math.abs(a.x - c[0]) > h[0] || Math.abs(a.y - c[1]) > h[1] || Math.abs(a.z - c[2]) > h[2]) return;
    if (DISP_BACKBONE.includes(a.atom) && a.resn !== 'GLY') return;
    (seen[a.chain] = seen[a.chain] || {})[a.resi] = a.resn;
  });
  return Object.keys(seen).map(ch => ({ chain: ch, resi: Object.keys(seen[ch]).map(Number), names: seen[ch] }));
}

function dispAddSticks(sel) {
  const p = state.proteinModel;
  if (DISP.colorBy === 'residue') {
    DISP_CLASSES.forEach(c => p.setStyle(Object.assign({}, sel, { resn: c.res }), { stick: { color: c.color, radius: 0.2 } }, true));
  } else p.setStyle(sel, { stick: Object.assign({ radius: 0.2 }, dispStickSpec()) }, true);
}

function dispSideChains() {
  if (DISP.repr === 'sticks' || DISP.repr === 'spheres') return;          // every atom is already drawn
  let mode = DISP.sideChains;
  if (mode === 'auto') mode = state.ligModel ? 'ligand' : 'pocket';
  if (mode === 'none') return;
  const p = state.proteinModel;
  if (mode === 'all') { dispAddSticks({}); return; }
  if (mode === 'ligand') {
    const lm = state.ligModel || DISP.previewModel;
    if (lm && typeof pocketSerials === 'function') {
      const ser = pocketSerials(p, lm, 5.0);
      if (ser && ser.length) { dispAddSticks({ serial: ser }); return; }
    }
  }
  dispPocketGroups().forEach(g => dispAddSticks({ chain: g.chain, resi: g.resi }));
}

function dispClearLabels() {
  const v = state.viewer; if (!v) return;
  DISP.labelObjs.forEach(l => { try { v.removeLabel(l); } catch (e) {} });
  DISP.labelObjs = [];
}
function dispLabels() {
  if (!DISP.labels || !DISP.showProtein) return;
  const v = state.viewer, p = state.proteinModel, want = {};
  dispPocketGroups().forEach(g => g.resi.forEach(r => { want[g.chain + ':' + r] = g.names[r]; }));
  let n = 0;
  p.selectedAtoms({ atom: 'CA' }).forEach(a => {
    if (!want[a.chain + ':' + a.resi] || n >= 40) return;
    n++;
    DISP.labelObjs.push(v.addLabel(`${a.resn}${a.resi}`, { position: { x: a.x, y: a.y, z: a.z }, fontSize: 11, fontColor: 'white',
      backgroundColor: dispResColor(a.resn) || '#8d99ae', backgroundOpacity: 0.75, showBackground: true, inFront: true }));
  });
}

// surface: drawn only for the pocket by default (a whole-protein surface can exhaust the GPU and blank the view)
function dispSurfaceStyle() {
  const opacity = Number(DISP.opacity);
  return {
    opacity: Number.isFinite(opacity) ? Math.max(0.1, Math.min(0.9, opacity)) : 0.5,
    transparent: true,
    depthWrite: false,
    colorscheme: { prop: 'resn', map: DISP_SURFACE_RESIDUE_COLORS },
  };
}
function dispSurface(tok) {
  const v = state.viewer, p = state.proteinModel;
  if (DISP.surface === 'off' || !DISP.showProtein) return;
  const style = dispSurfaceStyle();
  let sels = [];
  if (DISP.surface === 'whole') {
    if (p.selectedAtoms({}).length > 9000) {
      showToast('This protein is too big for a whole-protein surface — showing the pocket surface instead.', 'warning', '⚠️');
      DISP.surface = 'pocket'; dispSync();
    } else sels = [{}];
  }
  if (DISP.surface === 'pocket') sels = dispPocketGroups().map(g => ({ chain: g.chain, resi: g.resi }));
  sels.forEach(sel => {
    try {
      const sp = v.addSurface($3Dmol.SurfaceType.VDW, style, sel, sel);
      if (sp && typeof sp.then === 'function') {
        sp.then(id => {
          if (tok !== DISP.tok) { try { v.removeSurface(id); } catch (e) {} }
          else {
            DISP.surfaceIds.push(id);
            try { v.setSurfaceMaterialStyle(id, dispSurfaceStyle()); } catch (e) {}
            v.render();
          }
        })
          .catch(err => console.warn('Surface failed:', err));
      }
    } catch (e) {
      console.warn('Surface failed:', e);
      showToast('The surface could not be drawn — try "Pocket only".', 'warning', '⚠️');
    }
  });
}

function dispLigandStyles() {
  const lm = state.ligModel;
  if (lm) {
    if (typeof styleLigand === 'function') styleLigand(lm, state.sticksThick === undefined ? true : state.sticksThick);
    if (!DISP.showLigand) lm.setStyle({}, {});
    else if (!DISP.hydrogens) lm.setStyle({ elem: 'H' }, {});
  }
}

// ── the main renderer: every control calls this ──────────────────────────────
function renderProtein() {
  const v = state.viewer, p = state.proteinModel;
  if (!v || !p) return;
  const cv = document.querySelector('#viewer-container canvas');
  if (cv && !cv._dispHooked) {
    cv._dispHooked = true;
    cv.addEventListener('webglcontextlost', e => {
      e.preventDefault(); DISP.surface = 'off'; dispSync();
      showToast('The 3D view ran out of graphics memory — surface switched off. Try "Pocket only".', 'warning', '⚠️');
    });
    cv.addEventListener('webglcontextrestored', () => dispReload());
  }
  const tok = ++DISP.tok;
  dispClearLabels();
  try { v.removeAllSurfaces(); } catch (e) {}
  DISP.surfaceIds = [];
  try {
    if (!DISP.showProtein) p.setStyle({}, {});
    else {
      if (DISP.colorBy === 'residue') {
        p.setStyle({}, dispStyleObj({ color: '#8d99ae' }));
        DISP_CLASSES.forEach(c => p.setStyle({ resn: c.res }, dispStyleObj({ color: c.color })));
      } else p.setStyle({}, dispStyleObj(dispBaseSpec()));
      p.setStyle({ hetflag: true }, { stick: Object.assign({ radius: 0.18 }, DISP.colorBy === 'element' ? { colorscheme: 'Jmol' } : { colorscheme: 'yellowCarbon' }) });
      dispSideChains();
    }
    if (!DISP.hydrogens) p.setStyle({ elem: 'H' }, {});
    dispLabels();
    dispSurface(tok);
    dispLigandStyles();
  } catch (e) { console.warn('renderProtein:', e); }
  dispPreviewDraw(false);
  if (DISP.preview && !state.ligModel && (!DISP.previewData || DISP.previewData.key !== dispPreviewKey())) dispPreviewFetch();
  dispLegend();
  try { v.render(); } catch (e) {}
}
function styleProtein() { renderProtein(); }   // replaces the old function so every part of the app uses the new renderer

// ── ligand preview before docking ────────────────────────────────────────────
function dispPreviewKey() {
  const t = (state.targets && state.targets[state.selectedTarget]) || {};
  return [state.selectedTarget, (typeof getSMILES === 'function' ? getSMILES() : ''), (t.center || []).join(',')].join('|');
}
async function dispPreviewFetch() {
  const smiles = typeof getSMILES === 'function' ? getSMILES() : '', tid = state.selectedTarget;
  if (!DISP.preview || !smiles || !tid) { DISP.previewData = null; dispPreviewDraw(false); dispLegend(); return; }
  const key = dispPreviewKey();
  if (DISP.previewPending === key) return;
  DISP.previewPending = key;
  try {
    const d = await apiFetch('/api/insight/preview', { method: 'POST', body: JSON.stringify({ target_id: tid, smiles }) });
    if (dispPreviewKey() !== key) return;
    DISP.previewData = Object.assign({ key }, d);
  } catch (e) { DISP.previewData = { key, error: e.message }; }
  finally { if (DISP.previewPending === key) DISP.previewPending = ''; }
  dispPreviewDraw(true); dispLegend();
}
function dispPreviewDraw(zoom) {
  const v = state.viewer; if (!v) return;
  if (DISP.previewModel) { try { v.removeModel(DISP.previewModel); } catch (e) {} DISP.previewModel = null; }
  const d = DISP.previewData;
  if (!DISP.preview || state.ligModel || !state.proteinModel || !d || d.error || d.key !== dispPreviewKey()) { try { v.render(); } catch (e) {} return; }
  const m = v.addModel(d.sdf, 'sdf');
  DISP.previewModel = m;
  m.setStyle({}, { stick: { colorscheme: 'magentaCarbon', radius: 0.22 }, sphere: { colorscheme: 'magentaCarbon', scale: 0.25 } });
  if (!DISP.showLigand) m.setStyle({}, {}); else if (!DISP.hydrogens) m.setStyle({ elem: 'H' }, {});
  if (zoom && DISP.zoomedFor !== state.selectedTarget) { v.zoomTo({ model: m.getID() }); v.zoom(0.5); DISP.zoomedFor = state.selectedTarget; }
  v.render();
}
function dispPreviewNow() {
  DISP.preview = true; DISP.zoomedFor = ''; DISP.previewData = null; dispSync();
  if (state.ligModel) { showStructure(state.selectedTarget); return; }
  dispPreviewFetch();
}
function dispOnSmilesChanged() {
  if (!DISP.preview) return;
  const smiles = getSMILES();
  if (state.ligModel && state.lastResult && smiles !== state.lastResult.smiles) { showStructure(state.selectedTarget); return; }  // drop the old pose so the new preview shows
  if (!state.ligModel) dispPreviewFetch();
}
function dispSchedule() { clearTimeout(DISP.timer); DISP.timer = setTimeout(dispOnSmilesChanged, 700); }

// ── hydrogens: load the protonated receptor file only when asked ─────────────
async function fetchReceptorPdb(targetId) {
  state.pdbCache = state.pdbCache || {};
  const wantH = DISP.hydrogens, key = targetId + (wantH ? '|H' : '');
  if (state.pdbCache[key]) return state.pdbCache[key];
  if (wantH) {
    try {
      const r = await fetch(`/api/receptor/${targetId}/pdb_h`);
      if (r.ok) {
        const txt = await r.text();
        if (txt.split('\n').some(l => /^(ATOM|HETATM)/.test(l) && l.slice(76, 78).trim() === 'H')) { state.pdbCache[key] = txt; return txt; }
      }
    } catch (e) {}
    showToast('No protonated (hydrogen) structure is available for this target — showing heavy atoms only.', 'warning', '⚠️');
  }
  if (state.pdbCache[targetId]) return state.pdbCache[targetId];
  const res = await fetch(`/api/receptor/${targetId}/pdb`);
  if (!res.ok) throw new Error(`Could not load receptor (${res.status})`);
  const txt = await res.text();
  state.pdbCache[targetId] = txt;
  return txt;
}
function dispReload() { if (state.selectedTarget) showStructure(state.selectedTarget, { pose: state.shownPose || undefined }); }

// ── legend ───────────────────────────────────────────────────────────────────
function dispLegend() {
  const el = document.getElementById('display-legend'); if (!el) return;
  const sw = (c, t, tip) => `<span class="dl-item" title="${dh(tip || '')}"><i style="background:${c}"></i>${t}</span>`;
  const rows = [];
  let prot = '';
  if (!DISP.showProtein) prot = '<span style="color:var(--text-3)">hidden</span>';
  else if (DISP.colorBy === 'residue') prot = DISP_CLASSES.map(c => sw(c.color, c.name, c.note)).join('');
  else if (DISP.colorBy === 'element') prot = DISP_ELEMENTS.slice(0, 4).map(e => sw(e[1], e[0])).join('') + '<span class="dl-note">(cartoon is grey; atoms are coloured on sticks)</span>';
  else if (DISP.colorBy === 'chain') prot = '<span class="dl-note">each protein chain has its own colour</span>';
  else if (DISP.colorBy === 'spectrum') prot = '<span class="dl-bar"></span><span class="dl-note">blue = start (N-terminus) → red = end (C-terminus)</span>';
  else prot = sw('#9fb0c3', 'Protein');
  rows.push(`<div><b>Protein</b> ${prot}</div>`);
  if (state.proteinModel && state.proteinModel.selectedAtoms({ hetflag: true }).length)
    rows.push(`<div>${sw('#ffd60a', 'Cofactor / hetero group (FAD, heme…)', 'Helper molecules bound inside the protein')}</div>`);
  const lig = [];
  if (state.ligModel) lig.push(sw('#2af598', 'Docked ligand (carbons)'));
  else if (DISP.previewModel) lig.push(sw('#e056fd', 'Preview of your molecule — NOT docked yet', 'Placed at the centre of the grid box, in an arbitrary orientation'));
  if (lig.length) lig.push(sw('#3050F8', 'N') + sw('#FF0D0D', 'O') + sw('#FFFF30', 'S'));
  if (lig.length) rows.push(`<div><b>Ligand</b> ${lig.join('')}</div>`);
  const pd = DISP.previewData;
  if (DISP.preview && !state.ligModel && pd) {
    if (pd.error) rows.push(`<div class="dl-note">Preview unavailable: ${dh(pd.error)}</div>`);
    else rows.push(`<div class="dl-note">Your molecule: ${pd.heavy_atoms} heavy atoms, about ${pd.length} Å long · the grid box is ${pd.box_min_edge} Å at its narrowest → ${pd.length <= pd.box_min_edge ? '✅ it fits' : '⚠️ it may not fit — enlarge the box'}.</div>`);
  }
  const other = [];
  if (DISP.showBox) other.push(sw('#00ffcc', 'Grid box', 'The region the docking search explores'));
  if (DISP.hydrogens) other.push(sw('#FFFFFF', 'Hydrogens'));
  if (DISP.surface !== 'off') other.push(`<span class="dl-note">semi-transparent surface (${DISP.surface === 'pocket' ? 'pocket only' : 'whole protein'})</span>`);
  if (typeof GUIDE !== 'undefined' && GUIDE.inter && !GUIDE.inter.error && GUIDE.linesOn) {
    Object.keys(GUIDE_INTER).filter(k => GUIDE.inter.counts[k]).forEach(k => other.push(sw(GUIDE_INTER[k].color, GUIDE_INTER[k].name + ' (dashed line)', GUIDE_INTER[k].tip)));
  }
  if (other.length) rows.push(`<div><b>Other</b> ${other.join('')}</div>`);
  if (state.ligModel) {
    const interactions = Object.entries(POSE_INTERACTION_TYPES).map(([key, meta]) => {
      const count = state.poseInteractionData?.counts?.[key];
      return `<span class="dl-item" style="color:${meta.color}" title="Estimated from the selected pose's geometry"><i class="dl-dash"></i>${meta.name}${Number.isFinite(Number(count)) ? ` (${count})` : ''}</span>`;
    }).join('');
    const status = state.poseInteractionsLoading ? '<span class="dl-note">calculating…</span>'
      : state.poseInteractionError ? `<span class="dl-note" title="${dh(state.poseInteractionError)}">analysis unavailable</span>`
      : state.interactionsEnabled ? '' : '<span class="dl-note">hidden</span>';
    rows.push(`<div><b>Estimated pose interactions</b> ${interactions}${status}<span class="dl-note">dashed lines show geometry-based ligand–protein contacts</span></div>`);
  }
  el.innerHTML = rows.join('');
}

// ── controls ─────────────────────────────────────────────────────────────────
function dispSync() {
  const set = (id, v) => { const e = document.getElementById(id); if (e) { if (e.type === 'checkbox') e.checked = !!v; else e.value = v; } };
  set('dp-repr', DISP.repr); set('dp-color', DISP.colorBy); set('dp-side', DISP.sideChains); set('dp-surface', DISP.surface);
  set('dp-opacity', Math.round(DISP.opacity * 100)); set('dp-h', DISP.hydrogens); set('dp-labels', DISP.labels);
  set('dp-protein', DISP.showProtein); set('dp-ligand', DISP.showLigand); set('dp-box', DISP.showBox); set('dp-preview', DISP.preview);
  set('dp-interactions', state.interactionsEnabled);
  const interactionToggle = document.getElementById('dp-interactions');
  if (interactionToggle) interactionToggle.disabled = !state.ligModel;
}
function dispApplyBox() {
  const v = state.viewer; if (!v) return;
  if (!DISP.showBox) { (state.boxShapes || []).forEach(s => { try { v.removeShape(s); } catch (e) {} }); state.boxShapes = []; }
  else { const t = state.targets[state.selectedTarget]; if (t) drawGridBox(v, t.center, t.size, !!state.ligModel); }
}
function dispSet(key, val) {
  DISP[key] = val; dispSync();
  if (key === 'hydrogens') { dispReload(); return; }
  if (key === 'opacity' && state.viewer && DISP.surface !== 'off') {
    const style = dispSurfaceStyle();
    (DISP.surfaceIds || []).forEach(id => {
      try { state.viewer.setSurfaceMaterialStyle(id, style); } catch (e) {}
    });
    state.viewer.render();
    return;
  }
  if (key === 'showBox') dispApplyBox();
  if (key === 'preview' && val && !state.ligModel) dispPreviewFetch();
  renderProtein();
  if (typeof GUIDE !== 'undefined') {
    GUIDE.pocketOn = !!(DISP.labels && DISP.sideChains !== 'none');
    if (typeof guideSyncPocketBtn === 'function') guideSyncPocketBtn();
  }
}
function dispReset() {
  Object.assign(DISP, { repr: 'cartoon', colorBy: 'residue', sideChains: 'auto', surface: 'off', opacity: 0.5, labels: false,
    showProtein: true, showLigand: true, showBox: true, preview: true });
  const h = DISP.hydrogens; DISP.hydrogens = false; dispSync();
  if (h) dispReload(); else { dispApplyBox(); renderProtein(); }
}

// old toolbar buttons keep working
function toggleVis(type) {
  if (type === 'cartoon') { DISP.repr = 'cartoon'; DISP.surface = 'off'; dispSet('showProtein', true); }
  else if (type === 'surface') dispSet('surface', DISP.surface === 'off' ? 'pocket' : 'off');
  else if (type === 'sticks' && state.ligModel) { state.sticksThick = !(state.sticksThick === undefined ? true : state.sticksThick); renderProtein(); }
}
// the Guide's "show pocket chemistry" button now uses the same renderer
function guideTogglePocket() {
  const on = !(DISP.labels && DISP.sideChains !== 'none');
  DISP.labels = on; DISP.sideChains = on ? 'pocket' : 'auto';
  if (on) DISP.colorBy = 'residue';
  dispSync(); renderProtein();
  if (typeof GUIDE !== 'undefined') { GUIDE.pocketOn = on; if (typeof guideSyncPocketBtn === 'function') guideSyncPocketBtn(); }
}

function dispBuildUI() {
  const box = document.getElementById('viewer-container');
  if (!box || document.getElementById('display-panel')) return;
  const st = document.createElement('style');
  st.textContent = `#viewer-legend{display:none!important}
.disp-panel{margin:8px 0;border:1px solid var(--border);border-radius:10px;padding:8px 12px;background:var(--bg-1);font-size:12.5px;color:var(--text-2)}
.disp-panel summary{cursor:pointer;font-weight:700;color:var(--text-1,#fff)}
.disp-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:8px 14px;margin-top:8px}
.disp-grid label{display:flex;flex-direction:column;gap:3px;font-size:11.5px;color:var(--text-3)}
.disp-grid select,.disp-grid input[type=range]{width:100%}
.disp-checks{display:flex;flex-wrap:wrap;gap:6px 16px;margin-top:10px}
.disp-checks label{display:flex;align-items:center;gap:5px;color:var(--text-2)}
.disp-legend{margin-top:8px;padding:8px 12px;border:1px solid var(--border);border-radius:10px;background:var(--bg-1);font-size:12px;color:var(--text-2);display:flex;flex-direction:column;gap:5px}
.dl-item{display:inline-flex;align-items:center;gap:5px;margin:0 12px 0 8px;white-space:nowrap}
.dl-item i{display:inline-block;width:11px;height:11px;border-radius:3px;border:1px solid rgba(255,255,255,.35)}
.dl-dash{display:inline-block;width:22px;height:0!important;border:0!important;border-top:3px dashed currentColor!important;border-radius:0!important;vertical-align:middle}
.dl-note{color:var(--text-3);margin-left:8px}
.dl-bar{display:inline-block;width:120px;height:9px;border-radius:5px;vertical-align:middle;margin-left:8px;background:linear-gradient(90deg,#3050F8,#2ec4b6,#ffd60a,#ff4d6d)}`;
  document.head.appendChild(st);
  const opt = (v, t) => `<option value="${v}">${t}</option>`;
  const panel = document.createElement('details');
  panel.id = 'display-panel'; panel.className = 'disp-panel'; panel.open = true;
  panel.innerHTML = `<summary>🎛 Display options — choose what you want to see</summary>
    <div class="disp-grid">
      <label>Protein style<select id="dp-repr" onchange="dispSet('repr',this.value)">${opt('cartoon', 'Cartoon (ribbons)')}${opt('sticks', 'Sticks — every atom')}${opt('lines', 'Thin lines — every atom')}${opt('spheres', 'Spheres (space-filling)')}</select></label>
      <label>Colour by<select id="dp-color" onchange="dispSet('colorBy',this.value)">${opt('residue', 'Residue type (acidic / basic / polar…)')}${opt('element', 'Element (C N O S H)')}${opt('chain', 'Chain')}${opt('spectrum', 'Rainbow, start → end')}${opt('single', 'Single colour')}</select></label>
      <label>Show side chains<select id="dp-side" onchange="dispSet('sideChains',this.value)">${opt('auto', 'Auto (pocket; near ligand after docking)')}${opt('none', 'None — backbone only')}${opt('pocket', 'Lining the grid box')}${opt('ligand', 'Within 5 Å of the ligand')}${opt('all', 'Every residue')}</select></label>
      <label>Surface<select id="dp-surface" onchange="dispSet('surface',this.value)">${opt('off', 'Off')}${opt('pocket', 'Pocket only (fast, recommended)')}${opt('whole', 'Whole protein (slow)')}</select></label>
      <label>Surface opacity<input id="dp-opacity" type="range" min="10" max="90" step="5" value="50" oninput="dispSet('opacity',this.value/100)"></label>
    </div>
    <div class="disp-checks">
      <label><input id="dp-protein" type="checkbox" checked onchange="dispSet('showProtein',this.checked)">Protein</label>
      <label><input id="dp-ligand" type="checkbox" checked onchange="dispSet('showLigand',this.checked)">Ligand</label>
      <label><input id="dp-box" type="checkbox" checked onchange="dispSet('showBox',this.checked)">Grid box</label>
      <label><input id="dp-h" type="checkbox" onchange="dispSet('hydrogens',this.checked)">Hydrogens</label>
      <label><input id="dp-labels" type="checkbox" onchange="dispSet('labels',this.checked)">Label pocket residues</label>
      <label><input id="dp-interactions" type="checkbox" checked onchange="togglePoseInteractions(this.checked)">Show pose interaction lines</label>
      <label><input id="dp-preview" type="checkbox" checked onchange="dispSet('preview',this.checked)">Preview my molecule before docking</label>
      <button class="btn btn-sm" onclick="dispPreviewNow()">👁 Show my molecule now</button>
      <button class="btn btn-sm" onclick="dispReset()">↺ Reset display</button>
    </div>`;
  box.insertAdjacentElement('beforebegin', panel);
  const lg = document.createElement('div');
  lg.id = 'display-legend'; lg.className = 'disp-legend';
  box.insertAdjacentElement('afterend', lg);
  dispSync(); dispLegend();
}

// ── hook into the rest of the app (runs after all other scripts have loaded) ──
(function dispHooks() {
  const wrap = (name, fn) => { const o = window[name]; if (typeof o === 'function') window[name] = function () { const r = o.apply(this, arguments); try { fn.apply(this, arguments); } catch (e) { console.warn(e); } return r; }; };
  if (typeof window.showStructure === 'function') {
    const ss = window.showStructure;
    window.showStructure = function (tid, opts) {
      state.shownPose = (opts && opts.pose) || null; DISP.zoomedFor = ''; DISP.previewModel = null;
      const r = ss.apply(this, arguments);
      if (typeof GUIDE !== 'undefined') GUIDE.pocketOn = !!(DISP.labels && DISP.sideChains !== 'none');
      return r;
    };
  }
  if (typeof window.drawGridBox === 'function') {
    const dg = window.drawGridBox;
    window.drawGridBox = function () { if (DISP.showBox === false) return; return dg.apply(this, arguments); };
  }
  wrap('setSMILES', dispSchedule);
  const ub = window.updateBox;
  if (typeof ub === 'function') window.updateBox = async function () {
    const r = await ub.apply(this, arguments);
    DISP.previewData = null; renderProtein(); return r;
  };
  window.addEventListener('load', () => {
    dispBuildUI();
    const inp = document.getElementById('smiles-input');
    if (inp) inp.addEventListener('input', dispSchedule);
  });
})();
