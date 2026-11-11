// ══ CELL 35 · static/js/targets.js ══
/* ─── targets.js — targets — target cards, selecting a target, SMILES input (single source of truth) ─── */


// ── Loaded PDB details and structure selection ─────────────────────────────────
function renderTargets() {
  const grid = document.getElementById('target-grid');
  if (!grid) return;
  const tids = Object.keys(state.targets);
  if (!tids.length) {
    grid.innerHTML = '<div class="loaded-pdb-empty">No PDB structure is loaded. Import a PDB ID or upload a structure to begin.</div>';
    if (typeof clearProteinViewer === 'function') clearProteinViewer();
    return;
  }
  renderTargetStructureDetails(state.selectedTarget || null);
}

function summarizePdbCoordinates(pdbText) {
  const chains = new Set(), residues = new Set();
  let atomCount = 0;
  for (const line of String(pdbText || '').split(/\r?\n/)) {
    if (!line.startsWith('ATOM  ') && !line.startsWith('HETATM')) continue;
    if (line.startsWith('ATOM  ')) {
      atomCount++;
      const chain = line.slice(21, 22).trim() || '—';
      chains.add(chain);
      residues.add(`${chain}:${line.slice(22, 27).trim()}:${line.slice(17, 20).trim()}`);
    }
  }
  return { atomCount, chains: [...chains].sort(), residueCount: residues.size };
}

async function renderTargetStructureDetails(targetId) {
  const panel = document.getElementById('target-grid');
  const target = state.targets?.[targetId];
  if (!panel) return;
  if (!targetId || !target) {
    panel.innerHTML = '<div class="loaded-pdb-empty">Import a PDB ID or upload a structure to view its details here.</div>';
    return;
  }
  const requestId = state.targetDetailsRequest = (state.targetDetailsRequest || 0) + 1;
  const pdbId = String(target.pdb_id || targetId).toUpperCase();
  panel.innerHTML = `<div class="loaded-pdb-loading"><span class="spinner"></span><span>Loading details for ${escapeHtml(pdbId)}…</span></div>`;

  const proteinPromise = (async () => {
    state.pdbCache = state.pdbCache || {};
    if (state.pdbCache[targetId]) return state.pdbCache[targetId];
    const response = await fetch(`/api/receptor/${encodeURIComponent(targetId)}/pdb`);
    if (!response.ok) throw new Error(`PDB coordinates are not prepared (${response.status})`);
    const text = await response.text();
    state.pdbCache[targetId] = text;
    return text;
  })();
  const infoPromise = (async () => {
    state.targetInfoCache = state.targetInfoCache || {};
    if (!/^[A-Z0-9]{4}$/.test(pdbId)) return null;
    if (!state.targetInfoCache[pdbId]) {
      state.targetInfoCache[pdbId] = apiFetch(`/api/pdb/info/${encodeURIComponent(pdbId)}`).catch(() => null);
    }
    return state.targetInfoCache[pdbId];
  })();
  const [proteinResult] = await Promise.allSettled([proteinPromise]);
  if (state.targetDetailsRequest !== requestId || state.selectedTarget !== targetId) return;

  const stats = proteinResult.status === 'fulfilled' ? summarizePdbCoordinates(proteinResult.value) : null;
  const proteinError = proteinResult.status === 'rejected' ? proteinResult.reason?.message : '';
  const renderSummary = info => {
    if (state.targetDetailsRequest !== requestId || state.selectedTarget !== targetId) return;
    const center = (target.center || []).map(value => Number(value).toFixed(1)).join(', ');
    const size = (target.size || []).map(value => Number(value).toFixed(1)).join(' × ');
    const ligands = (state.nativeLigands?.length ? state.nativeLigands : target.native_ligands) || [];
    const ligandLabels = ligands.length
      ? ligands.map(lig => `${lig.ccd_id || lig.id || 'component'}${lig.name ? ` · ${lig.name}` : ''}${lig.chain && lig.chain !== '—' ? ` · chain ${lig.chain}` : ''}${lig.residue && lig.residue !== '—' ? ` · residue ${lig.residue}${lig.insertion || ''}` : ''}`).join(' · ')
      : (info?.binding_ligands || []).map(lig => `${lig.ccd_id || lig.id || 'component'}${lig.name ? ` · ${lig.name}` : ''}`).join(' · ')
        || target.native_ligand || target.site_ligand || 'Not listed';
    const details = [
      ['PDB ID', pdbId], ['Structure', info?.title || target.name || 'Protein structure'],
      ['Organism', info?.organism || 'Not listed'], ['Method', info?.method || 'Not listed'],
      ['Resolution', info?.resolution != null ? `${Number(info.resolution).toFixed(2)} Å` : 'Not listed'],
      ['Protein chains', stats ? (stats.chains.join(', ') || '—') : 'Coordinates not prepared'],
      ['Residues', stats ? stats.residueCount.toLocaleString() : '—'],
      ['Protein atoms', stats ? stats.atomCount.toLocaleString() : '—'],
      ['Native components', ligandLabels],
      ['Grid center (Å)', center || '—'], ['Grid dimensions (Å)', size || '—'],
    ];
    panel.innerHTML = `
      <article class="loaded-pdb-card">
        <header class="loaded-pdb-header">
          <div><div class="loaded-pdb-eyebrow">Selected protein structure</div><h3>${escapeHtml(info?.title || target.name || pdbId)}</h3></div>
          <span class="loaded-pdb-status ${target.ready ? 'ready' : 'prep'}">${target.ready ? 'Ready to dock' : 'Structure details'}</span>
        </header>
        <dl class="loaded-pdb-facts">${details.map(([label, value]) => `<div><dt>${escapeHtml(label)}</dt><dd>${escapeHtml(value)}</dd></div>`).join('')}</dl>
        ${proteinError ? `<div class="loaded-pdb-footnote">${escapeHtml(proteinError)}</div>` : ''}
      </article>`;
  };
  renderSummary(null);
  infoPromise.then(info => { if (info) renderSummary(info); }).catch(() => {});
}

function selectTarget(tid) {
  if (!tid || !state.targets?.[tid]) {
    state.selectedTarget = null;
    state.lbCurrentTab = null;
    if (typeof clearProteinViewer === 'function') clearProteinViewer();
    renderTargetStructureDetails(null);
    return;
  }
  state.selectedTarget = tid;
  state.lbCurrentTab = tid;
  renderLeaderboardTabs();
  refreshLeaderboard();

  // Show this target's grid box in the inputs
  const t = state.targets[tid];
  if (t) fillBoxInputs(t.center || [0, 0, 0], t.size || [22, 22, 22]);
  renderTargetStructureDetails(tid);
  renderPdbDockingActions(tid);
  if (typeof loadLiterature === 'function') loadLiterature(tid);

  // Auto-preview the protein in 3D viewer when selected
  if (typeof previewGridBox === 'function') {
    previewGridBox(false);
  }
}

let nativeLigandRequest = 0;
let nativeLigandController = null;
async function renderPdbDockingActions(targetId) {
  const panel = document.getElementById('pdb-docking-actions');
  if (!panel) return;
  const target = state.targets && state.targets[targetId];
  panel.hidden = !targetId || !target;
  if (panel.hidden) return;

  const select = document.getElementById('native-ligand-select');
  const choice = document.getElementById('native-ligand-choice');
  const sourceSelect = document.getElementById('docking-ligand-source');
  const note = document.getElementById('native-ligand-note');
  let initialNativeApplied = false;
  if (state.loadedNativeLigand && state.designerSmilesBeforeNative) {
    setSMILES(state.designerSmilesBeforeNative);
    state.designerSmilesBeforeNative = '';
  }
  state.loadedNativeLigand = null;
  if (select) {
    select.onchange = () => {
      if (select.value !== '') loadNativeLigand(Number(select.value));
    };
  }
  if (sourceSelect) sourceSelect.onchange = () => {
    if (sourceSelect.value === 'native') {
      if (select && select.value !== '') loadNativeLigand(Number(select.value));
    } else {
      state.loadedNativeLigand = null;
      if (state.designerSmilesBeforeNative) {
        setSMILES(state.designerSmilesBeforeNative);
        state.designerSmilesBeforeNative = '';
      }
      if (state.selectedTarget && typeof showStructure === 'function') showStructure(state.selectedTarget);
    }
  };
  const seq = ++nativeLigandRequest;
  if (nativeLigandController) nativeLigandController.abort();
  nativeLigandController = null;
  state.nativeLigands = [];
  if (select) select.replaceChildren();
  if (choice) choice.hidden = true;
  const nativeOption = sourceSelect?.querySelector('option[value="native"]');
  if (nativeOption) nativeOption.disabled = true;
  if (sourceSelect) sourceSelect.value = 'chosen';
  const applyLigands = candidates => {
    const previousNative = state.loadedNativeLigand;
    const previousSource = sourceSelect?.value || 'chosen';
    const previousSelected = select && select.value !== '' ? state.nativeLigands?.[Number(select.value)] : null;
    const ligands = Array.isArray(candidates) ? candidates : [];
    state.nativeLigands = ligands;
    renderTargetStructureDetails(targetId);
    if (select) {
      select.innerHTML = `<option value="" selected disabled>Choose native component…</option>` +
        ligands.map((lig, i) => {
          const atoms = lig.atom_count == null ? '' : ` · ${escapeHtml(lig.atom_count)} atoms`;
          const location = `${lig.chain && lig.chain !== '—' ? ` · chain ${escapeHtml(lig.chain)}` : ''}${lig.residue && lig.residue !== '—' ? ` · residue ${escapeHtml(lig.residue)}${escapeHtml(lig.insertion || '')}` : ''}`;
          return `<option value="${i}">${escapeHtml(lig.ccd_id)} · ${escapeHtml(lig.name || lig.ccd_id)}${atoms}${location}</option>`;
        }).join('');
    }
    if (ligands.length) {
      if (choice) choice.hidden = false;
      if (nativeOption) nativeOption.disabled = false;
      const nativeKey = lig => `${String(lig?.ccd_id || '').toUpperCase()}|${lig?.chain || ''}|${lig?.residue || ''}|${lig?.insertion || ''}`;
      const previousKey = previousNative ? nativeKey(previousNative) : (previousSelected ? nativeKey(previousSelected) : '');
      let selectedIndex = previousKey ? ligands.findIndex(lig => nativeKey(lig) === previousKey) : -1;
      if (selectedIndex < 0 && previousNative) {
        selectedIndex = ligands.findIndex(lig => String(lig.ccd_id || '').toUpperCase() === String(previousNative.ccd_id || '').toUpperCase());
      }
      let shouldLoadNative = previousSource === 'native' && !!previousNative && selectedIndex >= 0;
      if (selectedIndex < 0 && target.native_ligand_removed && !initialNativeApplied) {
        const preferred = String(target.native_ligand || target.site_ligand || '').toUpperCase();
        selectedIndex = ligands.findIndex(lig => String(lig.ccd_id || '').toUpperCase() === preferred);
        if (selectedIndex >= 0) {
          shouldLoadNative = true;
          initialNativeApplied = true;
        }
      }
      if (selectedIndex < 0 && previousSelected) {
        selectedIndex = ligands.findIndex(lig => String(lig.ccd_id || '').toUpperCase() === String(previousSelected.ccd_id || '').toUpperCase());
      }
      if (select && selectedIndex >= 0) select.value = String(selectedIndex);
      if (shouldLoadNative) loadNativeLigand(selectedIndex);
      const selectedNative = sourceSelect?.value === 'native' && state.loadedNativeLigand;
      if (note) note.textContent = selectedNative
        ? `${selectedNative.ccd_id} is selected at chain ${selectedNative.chain || '—'}${selectedNative.residue && selectedNative.residue !== '—' ? `, residue ${selectedNative.residue}${selectedNative.insertion || ''}` : ''}. The grid is centred on this ligand; choose another component to switch binding sites.`
        : 'Choose a PDB component below to dock it in its co-crystal site; selecting it also centres the grid box on that ligand.';
    } else if (note) {
      const component = target.native_ligand || target.site_ligand;
      note.textContent = component
        ? `This target records native ligand ${component}, but its ligand structure was not returned.`
        : 'No co-crystallized ligand was returned for this PDB. You can still dock the designer ligand.';
    }
  };
  const embeddedLigands = Array.isArray(target.native_ligands) ? target.native_ligands : [];
  if (embeddedLigands.length) applyLigands(embeddedLigands);
  if (note) note.textContent = 'Checking this PDB for co-crystallized ligands…';
  const controller = new AbortController();
  nativeLigandController = controller;
  const timeout = setTimeout(() => controller.abort(), 20000);
  try {
    const response = await apiFetch(`/api/targets/${encodeURIComponent(targetId)}/native-ligands`, { signal: controller.signal });
    if (seq !== nativeLigandRequest || state.selectedTarget !== targetId) return;
    const serverLigands = Array.isArray(response.ligands) ? response.ligands : [];
    const serverComponents = new Set(serverLigands.map(lig => String(lig.ccd_id || '').toUpperCase()));
    const embeddedFallbacks = embeddedLigands.filter(lig => !serverComponents.has(String(lig.ccd_id || '').toUpperCase()));
    applyLigands([...serverLigands, ...embeddedFallbacks]);
  } catch (error) {
    if (seq !== nativeLigandRequest || state.selectedTarget !== targetId) return;
    const component = target.native_ligand || target.site_ligand;
    const label = component ? ` The PDB lists component ${component}, but its 3-D chemistry could not be loaded.` : '';
    const reason = error.name === 'AbortError' ? 'the lookup timed out' : error.message;
    if (note) note.textContent = embeddedLigands.length
      ? `Showing the ligand already loaded from this PDB. Additional co-crystal components could not be checked: ${reason}.`
      : `Native PDB ligand unavailable: ${reason}.${label} You can still dock the designer ligand.`;
  } finally {
    clearTimeout(timeout);
    if (nativeLigandController === controller) nativeLigandController = null;
  }
}

function loadNativeLigand(index) {
  const native = state.nativeLigands?.[index];
  const targetId = state.selectedTarget;
  const target = state.targets && state.targets[targetId];
  if (!native || !targetId || !target) return;
  if (!state.loadedNativeLigand && !state.designerSmilesBeforeNative) state.designerSmilesBeforeNative = getSMILES();
  state.loadedNativeLigand = native;
  const sourceSelect = document.getElementById('docking-ligand-source');
  if (sourceSelect) sourceSelect.value = 'native';
  setSMILES(native.smiles || '');

  // Start at the co-crystal binding site, while leaving the coordinates and
  // box dimensions editable for a custom search.
  const inputs = readBoxInputs();
  const center = Array.isArray(native.center) && native.center.length === 3 ? native.center : inputs.center;
  const size = Array.isArray(native.size) && native.size.length === 3 ? native.size : inputs.size;
  target.center = center.slice(); target.size = size.slice();
  fillBoxInputs(center, size);
  renderTargetStructureDetails(targetId);
  updateBox();
  if (native.pdb && typeof showStructure === 'function') showStructure(targetId, { pose: native.pdb });
  else if (typeof previewGridBox === 'function') previewGridBox(false);
}

// ── SMILES input helpers ──────────────────────────────────────────────────────
// The #smiles-input text box is the single source of truth. The JSME editor
// writes into it whenever the drawing changes, and every loader (quick-fill,
// ligand browser, PubChem search, upload) writes into it too. Docking reads
// only the text box, so a loaded ligand can never be "lost" to a stale editor.
function getSMILES() {
  return document.getElementById('smiles-input').value.trim();
}

function setSMILES(smiles) {
  const clean = (smiles || '').trim();
  document.getElementById('smiles-input').value = clean;
  if (clean && window.jsmeApp && typeof window.jsmeApp.readSMILES === 'function') {
    state.jsmeSyncLock = true;              // don't let the editor echo overwrite us
    try { window.jsmeApp.readSMILES(clean); } catch (e) {}
    setTimeout(() => { state.jsmeSyncLock = false; }, 300);
  }
}
