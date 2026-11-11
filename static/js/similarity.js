/* ─── similarity.js — Drug Twin Finder: auto-discovery + manual compare (MACCS + Tanimoto made human) ─── */

const SIM_KNOWN = {
  MAOA:   [['Clorgyline (MAO-A drug)', 'C#CCN(C)CCCOc1ccc(Cl)cc1Cl'], ['Moclobemide (MAO-A drug)', 'O=C(NCCN1CCOCC1)c1ccc(Cl)cc1'], ['Harmine', 'COc1ccc2c(c1)[nH]c1c(C)nccc12']],
  MAOB:   [['Selegiline (MAO-B drug)', 'CC(Cc1ccccc1)N(C)CC#C'], ['Rasagiline (MAO-B drug)', 'C#CCNC1CCc2ccccc21'], ['Safinamide (MAO-B drug)', 'CC(NCc1ccc(OCc2cccc(F)c2)cc1)C(N)=O']],
  CYP3A4: [['Ketoconazole', 'CC(=O)N1CCN(CC1)c1ccc(OCC2COC(Cn3ccnc3)(O2)c2ccc(Cl)cc2Cl)cc1'], ['Midazolam', 'Cc1ncc2CN=C(c3ccccc3F)c3cc(Cl)ccc3-n12'],
           ['Nifedipine', 'COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC'], ['Testosterone', 'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O']],
};
const SIM = { timer: null, seq: 0, lastKey: '' };

function simPercentDisplay(row) {
  const percent = Number(row && row.percent);
  if (!Number.isFinite(percent)) return '—';
  if ((percent >= 100 && !row.fingerprint_identical) || (percent >= 99.995 && percent < 100)) return '>99.99%';
  return `${percent.toFixed(2)}%`;
}
function simPercentWidth(row) {
  const percent = Number(row && row.percent);
  return Number.isFinite(percent) ? Math.max(0, Math.min(100, percent)) : 0;
}

// ── Lipinski gates (Rule of Five) ────────────────────────────────────────────
function simLipBody() {
  const el = id => document.getElementById(id);
  return { lipinski1_enabled: el('sim-lip1-on') ? el('sim-lip1-on').checked : true,
           lipinski1_mode: el('sim-lip1-mode') ? el('sim-lip1-mode').value : 'standard',
           lipinski2_enabled: el('sim-lip2-on') ? el('sim-lip2-on').checked : true,
           lipinski2_mode: el('sim-lip2-mode') ? el('sim-lip2-mode').value : 'strict' };
}
function simLipChanged() { SIM.lastKey = ''; simSchedule(true); }

// "Will it work as a pill?" strip shown in the Chemistry Guide and on every twin card
function guideRo5HTML(f) {
  const l = f && (f.lipinski || f.ro5); if (!l) return '';
  const p = l.properties, c = l.checks || {};
  const item = (ok, label, val, lim) => `<span class="gd-chip" style="border-color:${ok ? '#2af598' : '#ff4d6d'};color:${ok ? '#2af598' : '#ff4d6d'}" title="Rule of Five limit: ${lim}">${ok ? '✅' : '❌'} ${label} ${val}</span>`;
  return `<div style="margin-top:8px"><b>💊 Lipinski's Rule of Five</b> — will it work as a pill?
    <span class="text-xs" style="color:var(--text-3)">(rule of thumb: more than one miss usually means poor oral absorption; many real drugs break it)</span>
    <div class="gd-chips" style="margin-top:4px">${item(c.mw, 'Weight', p.mw + ' Da', '≤ 500 Da')}${item(c.logp, 'logP', p.logp, '≤ 5')}${item(c.hbd, 'Donors', p.hbd, '≤ 5')}${item(c.hba, 'Acceptors', p.hba, '≤ 10')}</div>
    <div class="text-xs" style="margin-top:4px;color:var(--text-2)">${l.violations === 0 ? 'Passes all four limits.' : `${l.violations} miss${l.violations === 1 ? '' : 'es'}: ${l.violations_detail.map(sh).join('; ')}`}</div></div>`;
}
function simRo5Chip(r) {
  const l = r.ro5; if (!l) return '';
  const p = l.properties, ok = l.violations === 0, one = l.violations === 1;
  const col = ok ? '#2af598' : (one ? '#ffd60a' : '#ff4d6d');
  return `<span class="gd-chip" style="border-color:${col};color:${col}" title="${sh(l.violations ? l.violations_detail.join('; ') : 'Passes all four Rule-of-Five limits')}">💊 Ro5 ${ok ? '✅' : (one ? '⚠️' : '❌')} ${l.violations} miss${l.violations === 1 ? '' : 'es'} · MW ${Math.round(p.mw)} · logP ${p.logp} · HBD ${p.hbd} · HBA ${p.hba}</span>`;
}
function simGateSummary(data) {
  const g = data.lipinski; if (!g) return '';
  const n1 = g.gate1.rejected.length, n2 = g.gate2.rejected.length;
  const all = g.gate1.rejected.map(x => ({ ...x, gate: 1 })).concat(g.gate2.rejected.map(x => ({ ...x, gate: 2 })));
  let h = `<div class="text-xs" style="color:var(--text-2);margin:4px 0">💊 Lipinski gate 1 (${g.gate1.enabled ? g.gate1.mode : 'off'}): ${g.gate1.enabled ? `removed <b>${n1}</b> before comparing` : 'not applied'} · gate 2 (${g.gate2.enabled ? g.gate2.mode : 'off'}): ${g.gate2.enabled ? `removed <b>${n2}</b> after ranking` : 'not applied'}</div>`;
  if (all.length) h += `<details class="text-xs" style="color:var(--text-3);margin-bottom:6px"><summary>Why were ${all.length} molecule${all.length === 1 ? '' : 's'} filtered out?</summary>` +
    all.slice(0, 30).map(x => `<div>• <b>${sh(x.name)}</b> (gate ${x.gate}${x.percent ? ', ' + simPercentDisplay(x) + ' alike' : ''}): ${x.reasons.map(sh).join('; ')}</div>`).join('') +
    (all.length > 30 ? `<div>… and ${all.length - 30} more</div>` : '') + '<div>Untick a gate above to see them anyway.</div></details>';
  return h;
}
const sh = s => escapeHtml(String(s === undefined || s === null ? '' : s));

function simParseCandidates(text) {          // "Name, SMILES"  |  "SMILES Name"  |  "SMILES"
  const out = [];
  (text || '').split('\n').forEach(raw => {
    const line = raw.trim(); if (!line || line.startsWith('#')) return;
    let name = '', smiles = '';
    if (line.includes(',')) { const k = line.lastIndexOf(','); name = line.slice(0, k).trim(); smiles = line.slice(k + 1).trim(); }
    else { const parts = line.split(/\s+/); smiles = parts[0]; name = parts.slice(1).join(' '); }
    if (smiles) out.push({ name: name || `Compound ${out.length + 1}`, smiles });
  });
  return out;
}

function simUseDesigner() {
  const s = getSMILES();
  if (!s) { showToast('The designer is empty — draw or load a molecule first', 'warning', '⚠️'); return; }
  document.getElementById('sim-query').value = s;
}
function simAppend(list) {
  const box = document.getElementById('sim-candidates');
  const have = new Set(simParseCandidates(box.value).map(c => c.smiles));
  const add = list.filter(c => c.smiles && !have.has(c.smiles));
  if (!add.length) { showToast('Nothing new to add', 'info', 'ℹ️'); return; }
  box.value = (box.value.trim() ? box.value.trim() + '\n' : '') + add.map(c => `${c.name}, ${c.smiles}`).join('\n');
  showToast(`Added ${add.length} molecule(s)`, 'success', '🧪');
}
function simAddQuick() { simAppend(typeof QUICK_MOLS !== 'undefined' ? QUICK_MOLS : []); }
function simAddKnown() {
  const list = SIM_KNOWN[state.selectedTarget];
  if (!list) { showToast('No built-in reference drugs for this target — paste your own below (Name, SMILES)', 'info', 'ℹ️'); return; }
  simAppend(list.map(x => ({ name: x[0], smiles: x[1] })));
}
function simAddLastSearch() {
  const r = (state.lastLigandResults || []).filter(x => x.smiles);
  if (!r.length) { showToast('Run a PubChem search first', 'warning', '⚠️'); return; }
  simAppend(r.map(x => ({ name: x.name, smiles: x.smiles })));
}

// ── manual compare (the original "search" feature) ───────────────────────────
async function simRun() {
  const query = document.getElementById('sim-query').value.trim();
  const candidates = simParseCandidates(document.getElementById('sim-candidates').value);
  const box = document.getElementById('sim-results');
  if (!query) { showToast('Step 1: enter your molecule (or press "Use my designer molecule")', 'warning', '⚠️'); return; }
  if (!candidates.length) { showToast('Step 2: add at least one molecule to compare with', 'warning', '⚠️'); return; }
  box.innerHTML = '<div class="browser-placeholder"><div class="spinner" style="width:32px;height:32px"></div></div>';
  try {
    const data = await apiFetch('/api/similarity/compare', { method: 'POST', body: JSON.stringify(Object.assign({ query, candidates, fingerprint: 'maccs' }, simLipBody())) });
    box.innerHTML = '';
    simRender(data, box);
  } catch (err) { box.innerHTML = `<div style="color:var(--red);font-size:13px">Could not compare: ${sh(err.message)}</div>`; }
}

// ── AUTO-DISCOVERY: everything above the threshold, found for you ────────────
function simThreshold() { return (parseInt((document.getElementById('sim-th') || {}).value, 10) || 70) / 100; }

function simLocalLibrary() {
  const out = [], seen = new Set();
  const add = (name, smiles, src) => { if (smiles && !seen.has(smiles)) { seen.add(smiles); out.push({ name, smiles, source: src || 'library' }); } };
  (typeof QUICK_MOLS !== 'undefined' ? QUICK_MOLS : []).forEach(m => add(m.name, m.smiles));
  Object.keys(SIM_KNOWN).forEach(k => SIM_KNOWN[k].forEach(x => add(x[0], x[1])));
  (state.lastLigandResults || []).forEach(r => add(r.name, r.smiles, 'pubchem'));
  const ta = document.getElementById('sim-candidates');
  if (ta) simParseCandidates(ta.value).forEach(c => add(c.name, c.smiles, 'mine'));
  return out;
}

function simThChanged() {
  const v = document.getElementById('sim-th-val'); if (v) v.textContent = Math.round(simThreshold() * 100) + '%';
  simSchedule(true);
}
function simSchedule(force) {
  clearTimeout(SIM.timer);
  const on = document.getElementById('sim-auto-on');
  if (!force && on && !on.checked) return;
  SIM.timer = setTimeout(() => simDiscover(true), 1500);
}

// BEGIN REQUESTED LIGAND BROWSER SIMILARITY
// Dedicated ChEMBL search used by the Ligand Browser card. Keep these on
// window because the card calls them from inline HTML handlers.
window.updateChemblLigandThreshold = function () {
  const slider = document.getElementById('chembl-ligand-th');
  const label = document.getElementById('chembl-ligand-th-value');
  if (slider && label) label.textContent = `${slider.value}%`;
};

window.findSimilarDrugs = async function () {
  const status = document.getElementById('chembl-ligand-status');
  const box = document.getElementById('chembl-ligand-results');
  const button = document.querySelector('#chembl-ligand-browser-card button[onclick*="findSimilarDrugs"]');
  if (!status || !box) return;

  const query = String(typeof getSMILES === 'function' ? (getSMILES() || '') : '').trim();
  if (!query) {
    status.textContent = 'Load or draw a ligand in the Dock tab before searching ChEMBL.';
    box.innerHTML = '';
    return;
  }

  const slider = document.getElementById('chembl-ligand-th');
  const threshold = Math.max(0.4, Math.min(0.95, (Number(slider && slider.value) || 70) / 100));
  const buttonText = button ? button.textContent : '';
  if (button) { button.disabled = true; button.textContent = 'Searching ChEMBL…'; }
  status.textContent = `Searching ChEMBL for molecules at least ${Math.round(threshold * 100)}% similar…`;
  box.innerHTML = '<div class="browser-placeholder"><div class="spinner" style="width:24px;height:24px"></div></div>';

  try {
    const data = await apiFetch('/api/similarity/discover', {
      method: 'POST',
      body: JSON.stringify({
        query,
        threshold,
        use_chembl: true,
        use_pubchem: false,
        use_reference: false,
        local: [],
        max_results: 40,
        // A similarity lookup should still run for large or unusual ligands.
        lipinski1_enabled: false,
        lipinski1_mode: 'off',
        lipinski2_enabled: false,
        lipinski2_mode: 'off'
      })
    });
    const results = Array.isArray(data.results) ? data.results : [];
    const chemblNote = (data.notes || []).find(note => String(note).includes('ChEMBL:'));
    const count = Number(data.n_chembl) || 0;
    status.textContent = results.length
      ? `Found ${results.length} similar molecule${results.length === 1 ? '' : 's'}${count ? ` from ${count} ChEMBL compounds checked` : ''}.`
      : (chemblNote && String(chemblNote).includes('⚠️')
        ? `ChEMBL search failed: ${String(chemblNote).replace(/^.*ChEMBL:\s*/, '')}`
        : `No ChEMBL molecules met the ${Math.round(threshold * 100)}% similarity threshold${count ? ` (${count} ChEMBL compounds checked)` : ''}.`);
    box.innerHTML = '';
    if (results.length) {
      if (typeof window.simRender === 'function') window.simRender(data, box);
      else box.innerHTML = '<div class="text-xs" style="color:var(--red)">Search completed, but the result display did not load. Refresh the page and try again.</div>';
    }
  } catch (err) {
    status.textContent = `ChEMBL search failed: ${err && err.message ? err.message : String(err)}`;
    box.innerHTML = '';
  } finally {
    if (button) { button.disabled = false; button.textContent = buttonText || 'Find similar drugs in ChEMBL'; }
  }
};
// END REQUESTED LIGAND BROWSER SIMILARITY

async function simDiscover(auto) {
  const status = document.getElementById('sim-auto-status'), box = document.getElementById('sim-auto-results');
  if (!status || !box) return;
  const query = (typeof getSMILES === 'function' ? getSMILES() : '').trim();
  if (!query) { status.textContent = 'Draw or load a molecule in the designer and I will find its look-alikes automatically.'; box.innerHTML = ''; return; }
  const th = simThreshold(), useLib = document.getElementById('sim-src-lib').checked, usePc = document.getElementById('sim-src-pc').checked;
  const key = [query, th, useLib, usePc, JSON.stringify(simLipBody())].join('|');
  if (auto === true && key === SIM.lastKey) return;
  const seq = ++SIM.seq;
  status.innerHTML = `<span class="spinner" style="width:12px;height:12px;display:inline-block;vertical-align:middle"></span> Searching ${usePc ? 'PubChem and ' : ''}your libraries for molecules at least ${Math.round(th * 100)}% alike…`;
  try {
    const data = await apiFetch('/api/similarity/discover', { method: 'POST', body: JSON.stringify(Object.assign({
      query, threshold: th, use_pubchem: usePc, local: useLib ? simLocalLibrary() : [], max_results: 40, use_reference: useLib, use_chembl: usePc }, simLipBody())) });
    if (seq !== SIM.seq) return;
    SIM.lastKey = key;
    const n = data.results.length;
    status.innerHTML = `<b>${n}</b> look-alike${n === 1 ? '' : 's'} at ≥ ${Math.round(th * 100)}% · checked ${data.n_checked} molecules (${data.n_library} from your libraries, ${data.n_pubchem} from PubChem)` +
      (data.notes.length ? `<br><span style="color:var(--text-3)">${data.notes.map(sh).join('<br>')}</span>` : '') + simGateSummary(data);
    box.innerHTML = '';
    if (!n) box.innerHTML = `<div class="gd-sec">🆕 <b>Nothing at ${Math.round(th * 100)}% or above</b> — your molecule looks genuinely novel compared with everything checked. Lower the slider to see the nearest relatives.</div>`;
    else {
      const best = Math.max.apply(null, data.results.map(r => r.percent));
      box.insertAdjacentHTML('beforeend', `<div class="text-xs" style="color:var(--text-2);margin:6px 0">${best >= 90
        ? '⚠️ Almost the same structure as a known compound — probably not new.'
        : '💡 Related to known compounds — compare the 🟠 chips to see what is different about yours.'}</div>`);
      simRender(data, box);
    }
  } catch (err) {
    if (seq !== SIM.seq) return;
    status.innerHTML = `<span style="color:var(--red)">Could not search: ${sh(err.message)}</span>`;
  }
}

function simChips(list, kind) {
  const st = { shared: 'border-color:#2af598;color:#2af598', mine: 'border-color:#f8961e;color:#f8961e', theirs: 'border-color:#8d99ae;color:#8d99ae' }[kind];
  return list.map(f => `<span class="gd-chip" style="${st}">${f.emoji} ${sh(f.label)}</span>`).join('');
}

function simRender(data, box) {
  const lesson = r => r.percent >= 75 ? 'Structurally very close — similar behaviour is plausible, but similar structure never guarantees similar binding.'
    : (r.percent >= 55 ? 'Same family with some real differences — a reasonable starting point for ideas.'
    : (r.percent >= 35 ? 'Only loosely related, so the score alone tells you little.' : 'Low fingerprint overlap. That alone does not prove it is new, active or risky.'));
  const src = r => r.source === 'pubchem' ? (r.cid ? `<a href="https://pubchem.ncbi.nlm.nih.gov/compound/${r.cid}" target="_blank" rel="noopener" class="gd-chip" style="border-color:#4cc9f0;color:#4cc9f0;text-decoration:none">🌐 PubChem ${r.cid}</a>` : '<span class="gd-chip" style="border-color:#4cc9f0;color:#4cc9f0">🌐 PubChem</span>')
    : (r.source === 'mine' ? '<span class="gd-chip">✍️ your list</span>' : (r.source ? '<span class="gd-chip">📚 library</span>' : ''));
  let html = `<div class="text-xs" style="color:var(--text-2);margin-bottom:8px">Your molecule has <b>${data.query.maccs_bits}</b> of the 166 checklist boxes ticked. Closest relatives first:</div>`;
  if (box.id === 'sim-results') html += simGateSummary(data);
  data.results.forEach(r => {
    html += `<div class="gd-sec" style="margin-bottom:10px">
      <div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap">
        <b style="font-size:15px">#${r.rank} ${sh(r.name)}</b>
        <span class="gd-chip" style="border-color:#ffd60a;color:#ffd60a">${r.level_label}</span>
        <span style="color:#ffd60a">${'★'.repeat(r.stars)}${'☆'.repeat(5 - r.stars)}</span>${typeof simSourceBadge === 'function' ? simSourceBadge(r) : src(r)}
        <span style="margin-left:auto"><b style="font-size:20px">${simPercentDisplay(r)}</b> alike</span>
        <button class="btn btn-sm sim-load" data-i="${r.index}">⚡ Load to dock</button></div>
      <div style="height:8px;border-radius:4px;background:rgba(255,255,255,.08);margin:8px 0;overflow:hidden">
        <div style="width:${simPercentWidth(r)}%;height:100%;background:linear-gradient(90deg,#4cc9f0,#f72585)"></div></div>
      <div class="gd-chips">${simChips(r.shared_features, 'shared')}${simChips(r.only_query, 'mine')}${simChips(r.only_candidate, 'theirs')}</div>
      <div class="gd-chips" style="margin-top:6px">${simRo5Chip(r)}${typeof simSynthChip === 'function' ? simSynthChip(r) : ''}</div>
      <div class="text-xs" style="color:var(--text-3);margin:4px 0">🟢 you both have it · 🟠 only YOUR molecule has it · ⚪ only THEIRS has it</div>
      <div class="text-xs" style="color:var(--text-2)">💡 ${lesson(r)}</div>
      <details class="text-xs" style="color:var(--text-3);margin-top:4px"><summary>Nerd view</summary>
        MACCS Tanimoto = ${r.maccs} = shared boxes ÷ (yours + theirs − shared) = ${r.bits_shared} ÷ (${r.bits_query} + ${r.bits_candidate} − ${r.bits_shared}).
        Finer ECFP4 fingerprint: ${r.morgan} (it is always lower — it looks at atom neighbourhoods, not just feature types).</details></div>`;
  });
  if (!data.results.length) html = '<div style="color:var(--text-3)">No results.</div>';
  if (data.invalid && data.invalid.length) html += `<div class="text-xs" style="color:var(--red)">Skipped (invalid SMILES): ${data.invalid.map(sh).join(', ')}</div>`;
  box.insertAdjacentHTML('beforeend', html);
  box.querySelectorAll('.sim-load').forEach(btn => btn.addEventListener('click', () => {
    const r = data.results.find(x => String(x.index) === btn.dataset.i);
    if (r) loadLigand(r.smiles, r.name);
  }));
}

// badge on every PubChem search result: how alike is it to the molecule in the designer?
async function simBadgeResults(results) {
  const q = (typeof getSMILES === 'function') ? getSMILES() : '';
  if (!q || !results || !results.length) return;
  try {
    const cands = results.map(r => ({ name: r.name, smiles: r.smiles || '' })).filter(c => c.smiles);
    const data = await apiFetch('/api/similarity/compare', { method: 'POST', body: JSON.stringify(Object.assign({ query: q, candidates: cands, fingerprint: 'maccs', top_n: 500 }, simLipBody())) });
    const by = {}; data.results.forEach(r => { by[r.smiles] = r; });
    const cards = document.querySelectorAll('#ligand-results .browser-result-card');
    results.forEach((r, i) => {
      const hit = by[(r.smiles || '').trim()], meta = cards[i] && cards[i].querySelector('.result-meta');
      if (hit && meta) meta.insertAdjacentHTML('beforeend', `<span class="result-tag" title="${hit.same ? 'Exact structure match' : 'Fingerprint similarity; not an exact identity'}">${hit.level_label} ${simPercentDisplay(hit)}</span>`);
    });
  } catch (e) { console.warn('similarity badges skipped:', e); }
}

// ── auto-run whenever the designer molecule changes (or the card scrolls into view) ──
(function simHooks() {
  const o = window.setSMILES;
  if (typeof o === 'function') window.setSMILES = function () { const r = o.apply(this, arguments); try { simSchedule(false); } catch (e) {} return r; };
  window.addEventListener('load', () => {
    const inp = document.getElementById('smiles-input');
    if (inp) inp.addEventListener('input', () => simSchedule(false));
    const card = document.getElementById('similarity-card');
    if (card && typeof IntersectionObserver !== 'undefined')
      new IntersectionObserver(es => { if (es.some(e => e.isIntersecting)) { const on = document.getElementById('sim-auto-on'); if (on && on.checked) simSchedule(true); } }).observe(card);
  });
})();

// BEGIN MANUAL_SIMILARITY_V1
(function () {
  'use strict';

  // The existing “Find look-alikes now” button still works.
  // Only automatic searches after editing/loading are disabled.
  window.simSchedule = function () {
    const status = document.getElementById('sim-auto-status');
    if (status) status.textContent = 'Molecule changed. Press “Find look-alikes now” when you want to search.';
  };
})();
// END MANUAL_SIMILARITY_V1
