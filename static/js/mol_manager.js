/* mol_manager.js — Cell 49p: drawing library · universal auto-render into the 2D editor · pinned reference docking score */
(function () {
  'use strict';
  if (window.__mm49p) return;
  window.__mm49p = true;

  const $ = id => document.getElementById(id);
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const toast = (m, t, i) => { try { showToast(m, t || 'info', i || '🧪'); } catch (e) {} };
  const boxSmiles = () => { const i = $('smiles-input'); return i ? i.value.trim() : ''; };
  const sync = ms => { try { state.jsmeSyncLock = true; setTimeout(() => { state.jsmeSyncLock = false; }, ms || 600); } catch (e) {} };
  const jsme = () => (window.jsmeApp && (typeof window.jsmeApp.readMolFile === 'function' || typeof window.jsmeApp.readGenericMolecularInput === 'function')) ? window.jsmeApp : null;
  const jsmeSmi = () => { try { const j = jsme(); const s = j ? String(j.smiles() || '').trim() : ''; return s === '.' ? '' : s; } catch (e) { return ''; } };
  const short = s => (s.length > 60 ? s.slice(0, 57) + '…' : s);
  function wrap(name, mk) {
    const o = window[name];
    if (typeof o !== 'function' || o.__mm) return;
    const w = mk(o); w.__mm = true; window[name] = w;
  }

  // ══ 1) Universal auto-render into the 2D editor ═══════════════════════════════
  let drawSeq = 0;                                       // only the newest load may draw (fast clicking never shows a stale molecule)
  function putGeneric(j, smi) {                          // fallback when the server cannot be reached
    try { if (typeof j.readGenericMolecularInput === 'function') { sync(1200); j.readGenericMolecularInput(smi); return true; } } catch (e) {}
    try { if (typeof j.readSMILES === 'function') { sync(1200); j.readSMILES(smi); return true; } } catch (e) {}
    return false;
  }
  async function viaMolfile(smi, seq) {                  // RDKit makes the 2D drawing (and checks the SMILES is valid)
    try {
      const d = await apiFetch('/api/chem/molfile?smiles=' + encodeURIComponent(smi));
      if (seq !== undefined && seq !== drawSeq) return true;          // a newer molecule was loaded meanwhile
      const j = jsme();
      if (j && d && d.molfile) { sync(1200); j.readMolFile(d.molfile); return true; }
    } catch (e) {
      if (/not a valid smiles/i.test(String(e && e.message))) return 'invalid';   // half-typed SMILES: leave the editor alone
    }
    return false;
  }
  function renderMol(smi, molfile, tries) {              // waits for the editor to finish loading, then draws
    smi = (smi || '').trim(); tries = tries || 0;
    if (!smi && !molfile) return;
    const j = jsme();
    if (!j) { if (tries < 120) setTimeout(() => renderMol(smi, molfile, tries + 1), 250); return; }
    const seq = ++drawSeq;
    if (molfile) { sync(1200); try { j.readMolFile(molfile); } catch (e) { if (smi) putGeneric(j, smi); } return; }
    viaMolfile(smi, seq).then(ok => { if (ok === false && seq === drawSeq) { const jj = jsme(); if (jj) putGeneric(jj, smi); } });
  }
  function ensureShown(smi) { if (smi) renderMol(smi); }   // called after the app's own setSMILES
  window.renderMoleculeInJSME = (smi, molfile) => renderMol(smi, molfile);

  // ══ 2) Pinned reference ═══════════════════════════════════════════════════════
  // Pin = the molecule (and protein) the "score to beat" belongs to. Editing the drawing never touches it.
  let Pin = { smiles: '', name: '', tid: '', cache: null, next: null };
  window.MMPin = () => Pin;
  let keep = false;                                      // true while loading one of YOUR saved drawings (= a modification)

  function newReference(smi, name) {                     // a completely different molecule was loaded
    if (!smi) return;
    Pin.next = { smiles: smi, name: name || 'the loaded molecule' };
    Pin.cache = null;
  }

  async function pinnedRefresh(start, force) {           // replaces plRefresh from Cell 49o
    if (!PL.analysis) return;
    if (Pin.next || !Pin.smiles || Pin.tid !== PL.tid) { // (re)pin: first load, new molecule, new protein or reset button
      const nm = (Pin.next && Pin.next.name) || (Pin.smiles && Pin.smiles === PL.smiles ? Pin.name : '') || 'the loaded molecule';
      Pin = { smiles: PL.smiles, name: nm, tid: PL.tid, cache: null, next: null };
    }
    const mine = Pin, env = plEng() + '|' + plExh();
    const c = mine.cache;
    if (!force && c && c.env === env && c.data.best && !c.data.pending.length) {      // complete benchmark → frozen, no server call
      Object.assign(PL, c.data, { error: '' });
      plRenderRef(); plRenderBanner(); plEvaluate();
      return;
    }
    try {
      const d = await plPost('/api/pipeline/best_known', { smiles: mine.smiles, target_id: mine.tid, engine: plEng(), exhaustiveness: plExh(), player: plPlayer() });
      if (mine !== Pin || env !== plEng() + '|' + plExh()) return;                  // protein / reference / engine changed meanwhile
      const data = { rows: d.rows, best: d.best, pending: d.pending, bestEngine: d.engine, loose: d.loose, same: d.same, rationale: d.rationale };
      mine.cache = { env, data };
      Object.assign(PL, data, { error: '' });
    } catch (e) { if (mine !== Pin) return; PL.error = e.message; }
    plRenderRef(); plRenderBanner(); plEvaluate();
    if (start && plAuto()) plDockPending();              // docks whatever known scores are still missing
  }
  window.plRefresh = pinnedRefresh;

  window.plRedockAll = async function () {               // re-dock the benchmark molecules (not the edited molecule)
    if (!PL.analysis) return;
    const seq = ++PL.dockSeq;
    for (const r of (PL.rows || []).filter(x => x.kind !== 'community')) {
      if (seq !== PL.dockSeq) return;
      PL.docking = r.smiles; PL.error = ''; plRenderRef(); plRenderBanner();
      try { await plPost('/api/pipeline/reference_dock', { target_id: PL.tid, smiles: r.smiles, name: r.name, exhaustiveness: plExh(), engine: plEng(), refresh: true }); }
      catch (e) { PL.error = r.name + ': ' + e.message; break; }
      await pinnedRefresh(false, true);
    }
    if (seq === PL.dockSeq) { PL.docking = ''; plRenderRef(); plRenderBanner(); }
  };

  window.mmResetReference = function () {
    const smi = boxSmiles();
    if (!smi) { toast('The editor is empty — load or draw a molecule first.', 'warning', '⚠️'); return; }
    Pin.next = { smiles: smi, name: 'your molecule' }; Pin.cache = null;
    try { PL.dockSeq++; PL.docking = ''; } catch (e) {}
    toast('Reference reset: the score to beat is now found for the molecule in the editor.', 'info', '🔄');
    plAnalyze();
  };

  const origRenderRef = window.plRenderRef;
  window.plRenderRef = function () {
    origRenderRef.apply(this, arguments);
    try {
      const box = $('pl-ref'); if (!box || !PL.analysis || !Pin.smiles) return;
      const same = boxSmiles() === Pin.smiles;
      const d = document.createElement('div');
      d.id = 'mm-pin'; d.style.cssText = 'margin-top:8px;padding:8px;border:1px dashed var(--border,#4a5568);border-radius:6px';
      d.innerHTML = '<div class="text-xs"><b>🔒 Pinned reference:</b> ' + esc(Pin.name) + ' on <b>' + esc(plProtein()) +
        '</b>. Edit your molecule freely: this score to beat stays fixed. It resets when you change protein, load a different molecule, or press the button.</div>' +
        '<div class="text-xs" style="margin:4px 0;color:var(--text-3)">Reference SMILES: ' + esc(short(Pin.smiles)) +
        (same ? ' (same as the editor)' : ' · you are working on a modified molecule') + '</div>' +
        '<button class="btn btn-sm" onclick="mmResetReference()">🔄 Reset reference to the molecule in the editor</button>';
      box.appendChild(d);
    } catch (e) {}
  };

  // ══ 3) Hooks: every way of loading a molecule ═════════════════════════════════
  let explicitName = null;
  document.addEventListener('click', ev => {             // quick-fill chips
    const el = ev.target && ev.target.closest && ev.target.closest('.qf-chip');
    if (el) { explicitName = (el.textContent || '').trim() || 'quick-fill molecule'; setTimeout(() => { explicitName = null; }, 100); }
  }, true);

  wrap('setSMILES', o => function (smiles) {
    const clean = (smiles || '').trim();
    const named = explicitName; explicitName = null;
    if (!keep && clean && (named || clean !== boxSmiles())) newReference(clean, named);   // initJSME re-sending the same text is not a new drug
    const r = o.apply(this, arguments);
    ensureShown(clean);
    return r;
  });
  wrap('loadLigand', o => function (smiles, name) {      // PubChem browser / search / Drug Twin Finder
    explicitName = name || 'loaded molecule';
    try { return o.apply(this, arguments); } finally { explicitName = null; }
  });
  wrap('selectTarget', o => function (tid) {             // new protein → reference resets automatically
    if (Pin.tid && tid !== Pin.tid) {
      Pin.next = { smiles: boxSmiles(), name: boxSmiles() === Pin.smiles ? Pin.name : 'your molecule' }; Pin.cache = null;
      toast('Protein changed — the reference score was reset for the new protein.', 'info', '🔄');
    }
    return o.apply(this, arguments);
  });
  window.handleLigandUpload = function (ev) {            // .smi / .mol / .sdf → drawn in the editor, becomes the new reference
    const file = ev.target.files && ev.target.files[0]; if (!file) return;
    const base = file.name.replace(/\.[^/.]+$/, '');
    const rd = new FileReader();
    rd.onload = e => {
      const text = String(e.target.result || '');
      if (/\.(smi|smiles|txt)$/i.test(file.name)) {
        const first = (text.split(/\r?\n/).find(l => l.trim()) || '').trim().split(/\s+/)[0];
        if (!first) { toast('Could not read a SMILES from that file', 'error', '⚠️'); return; }
        explicitName = base; window.setSMILES(first);
      } else {
        const j = jsme(); if (!j) { toast('The 2D editor is not ready yet', 'error', '⚠️'); return; }
        sync(1000);
        try { j.readMolFile(text.split('$$$$')[0]); } catch (er) { toast('Could not parse that file', 'error', '⚠️'); return; }
        const s = jsmeSmi();
        if (!s) { toast('Could not parse ligand', 'error', '⚠️'); return; }
        explicitName = base; window.setSMILES(s);
      }
      toast('Loaded ' + file.name + ' into the drawing interface', 'success', '🧪');
    };
    rd.readAsText(file); ev.target.value = '';
  };

  function bindTyping() {                                // typed or pasted SMILES also appear in the editor
    const inp = $('smiles-input'); if (!inp || inp.__mm) return; inp.__mm = true;
    let t = null;
    inp.addEventListener('input', ev => {
      if (!ev.isTrusted) return;
      if (/paste|drop/i.test(ev.inputType || '')) newReference(boxSmiles(), 'Pasted molecule');   // a pasted SMILES = a different molecule
      clearTimeout(t);
      t = setTimeout(() => { const v = boxSmiles(); if (v) viaMolfile(v, ++drawSeq); }, 700);                // invalid / half-typed SMILES is ignored
    });
  }

  // ══ 4) Drawing library ════════════════════════════════════════════════════════
  const LS = 'medchem_saved_drawings';
  let drawings = [];
  const lsRead = () => { try { return JSON.parse(localStorage.getItem(LS) || '[]'); } catch (e) { return []; } };
  const lsWrite = () => { try { localStorage.setItem(LS, JSON.stringify(drawings)); } catch (e) {} };

  function refreshSelect(selId) {
    const sel = $('mm-drawings'); if (!sel) return;
    let h = '<option value="">📂 Load saved drawing (' + drawings.length + ')…</option>';
    drawings.forEach(d => { h += '<option value="' + esc(d.id) + '"' + (d.id === selId ? ' selected' : '') + '>' + esc(d.name) + ' · ' + esc(d.created || '') + '</option>'; });
    sel.innerHTML = h;
  }
  async function loadDrawings() {
    const local = lsRead();
    try {
      const d = await apiFetch('/api/drawings'); drawings = d.drawings || [];
      for (const it of local) {                          // move drawings saved only in this browser up to the server
        if (drawings.some(x => x.id === it.id) || !it.smiles) continue;
        try { await apiFetch('/api/drawings', { method: 'POST', body: JSON.stringify(it) }); drawings.push(it); } catch (e) {}
      }
    } catch (e) { drawings = local; }
    lsWrite(); refreshSelect();
  }

  window.mmSaveDrawing = async function () {
    let smi = jsmeSmi() || boxSmiles(), molfile = '';
    try { const j = jsme(); if (j) molfile = j.molFile() || ''; } catch (e) {}
    if (!smi) { toast('The editor is empty. Draw or load a molecule first!', 'warning', '⚠️'); return; }
    const name = (prompt('Name for this drawing:', 'Drawing ' + (drawings.length + 1)) || '').trim();
    if (!name) return;
    const created = new Date().toLocaleString([], { year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' });
    let item = { id: 'l_' + Date.now(), name, smiles: smi, molfile, created };
    try { item = await apiFetch('/api/drawings', { method: 'POST', body: JSON.stringify(item) }); }
    catch (e) { toast('Saved in this browser only (' + e.message + ')', 'warning', '💾'); }
    drawings = [item].concat(drawings.filter(d => d.id !== item.id)); lsWrite(); refreshSelect(item.id);
    toast('Saved "' + name + '" to your drawing library', 'success', '💾');
  };

  window.mmLoadDrawing = function (id) {
    if (!id) return;
    const f = drawings.find(d => d.id === id); if (!f) return;
    keep = true;                                         // your own saved molecule counts as a modification: the reference stays pinned
    try { window.setSMILES(f.smiles); } finally { keep = false; }
    renderMol(f.smiles, f.molfile || '');                // exact drawing (layout included)
    toast('Loaded "' + f.name + '" into the drawing interface', 'success', '📂');
  };

  window.mmDeleteDrawing = async function () {
    const sel = $('mm-drawings');
    if (!sel || !sel.value) { toast('Pick a saved drawing in the dropdown first.', 'info', 'ℹ️'); return; }
    const it = drawings.find(d => d.id === sel.value); if (!it) return;
    if (!confirm('Delete "' + it.name + '" from your saved drawings?')) return;
    try { await apiFetch('/api/drawings/' + encodeURIComponent(it.id), { method: 'DELETE' }); } catch (e) {}
    drawings = drawings.filter(d => d.id !== it.id); lsWrite(); refreshSelect();
    toast('Deleted "' + it.name + '"', 'info', '🗑️');
  };

  window.mmExportDrawing = function () {
    let molfile = '', smi = jsmeSmi() || boxSmiles();
    try { const j = jsme(); if (j) molfile = j.molFile() || ''; } catch (e) {}
    if (!smi) { toast('Nothing to export.', 'warning', '⚠️'); return; }
    const ext = molfile ? 'mol' : 'smi';
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([molfile || smi], { type: 'text/plain;charset=utf-8' }));
    a.download = 'molecule_' + Date.now() + '.' + ext; a.click();
    toast('Exported the drawing as .' + ext, 'success', '⬇️');
  };

  function buildToolbar() {
    if ($('mm-bar')) return;
    const wrapEl = document.querySelector('.jsme-wrap') || ($('jsme-container') && $('jsme-container').parentElement);
    if (!wrapEl) return;
    const bar = document.createElement('div');
    bar.id = 'mm-bar'; bar.className = 'flex items-center gap-2 mt-2';
    bar.style.cssText = 'flex-wrap:wrap;display:flex;gap:6px;align-items:center';
    bar.innerHTML = '<button class="btn btn-sm btn-primary" onclick="mmSaveDrawing()" title="Save what is drawn in the editor">💾 Save drawing</button>' +
      '<select id="mm-drawings" class="smiles-input" style="width:auto;min-width:170px;max-width:240px;padding:4px 6px;font-size:11px" onchange="mmLoadDrawing(this.value)" title="Load a saved drawing back into the editor"><option value="">📂 Load saved drawing…</option></select>' +
      '<button class="btn btn-sm" onclick="mmDeleteDrawing()" title="Delete the selected saved drawing" style="font-size:11px">🗑️</button>' +
      '<button class="btn btn-sm" onclick="mmExportDrawing()" title="Download the drawing as a .mol file" style="font-size:11px">⬇️ Export .mol</button>';
    wrapEl.insertAdjacentElement('afterend', bar);
  }

  function init() { buildToolbar(); bindTyping(); loadDrawings(); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
