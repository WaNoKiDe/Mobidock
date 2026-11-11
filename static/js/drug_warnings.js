/* Route-specific RDKit warning panel; no aggregate drug-likeness score. */
(function () {
  'use strict';
  const W = { timer: null, seq: 0, smiles: '', data: null, site: 'oral' };
  const esc = value => escapeHtml(String(value == null ? '' : value));
  const $ = id => document.getElementById(id);

  function render() {
    const box = $('drug-warning-output');
    if (!box || !W.data) return;
    const all = W.data, site = all.by_site[W.site] ? W.site : all.default;
    W.site = site;
    const w = all.by_site[site];
    const options = all.order.map(id => `<option value="${id}"${id === site ? ' selected' : ''}>${esc(all.labels[id])}</option>`).join('');
    const list = (w.groups || []).map(group => `<section style="margin-top:8px"><b>${esc(w.by_group[group][0].group_label)}</b>` +
      w.by_group[group].map(item => `<div class="text-xs" style="margin:4px 0 0 8px;color:${item.severe ? '#ff4d6d' : 'var(--text-2)'}">${item.severe ? '⛔' : '⚠️'} ${esc(item.text)}</div>`).join('') + '</section>').join('');
    const color = w.blocking.length ? '#ff4d6d' : (w.count ? '#ffd60a' : '#2af598');
    box.innerHTML = `<div class="gd-sec" style="margin-top:8px"><label class="text-xs" style="display:block">Intended site of action
      <select id="drug-warning-site" class="smiles-input" style="width:100%;margin-top:4px">${options}</select></label>
      <div style="font-weight:700;color:${color};margin-top:8px">${esc(w.verdict)}</div>
      <div class="text-xs" style="color:var(--text-3);margin-top:4px">${w.descriptors_checked} property limits checked for this use; ${all.rules} rules across ${all.order.length} intended uses.</div>
      ${list || '<div class="text-xs" style="margin-top:8px;color:var(--text-2)">No checked property threshold was crossed for this use.</div>'}
      <div class="text-xs" style="color:var(--text-3);margin-top:8px">${esc(all.note)}</div></div>`;
    $('drug-warning-site').addEventListener('change', event => { W.site = event.target.value; render(); });
  }

  async function update(force) {
    const box = $('drug-warning-output');
    const smiles = (typeof getSMILES === 'function' ? getSMILES() : '').trim();
    if (!box) return;
    if (!smiles) { W.smiles = ''; W.data = null; box.innerHTML = '<div class="text-xs" style="color:var(--text-3);padding:10px">Draw or load a ligand to see route-specific property warnings.</div>'; return; }
    if (!force && smiles === W.smiles && W.data) return;
    const seq = ++W.seq;
    box.innerHTML = '<div class="text-xs" style="padding:10px;color:var(--text-3)">Checking RDKit property thresholds…</div>';
    try {
      const data = await apiFetch('/api/drug-warnings', { method: 'POST', body: JSON.stringify({ smiles }) });
      if (seq !== W.seq) return;
      W.smiles = smiles; W.data = data; render();
    } catch (error) {
      if (seq === W.seq) box.innerHTML = `<div class="text-xs" style="padding:10px;color:var(--red)">Warnings unavailable: ${esc(error.message)}</div>`;
    }
  }
  function schedule(force) { clearTimeout(W.timer); W.timer = setTimeout(() => update(force === true), 300); }

  const original = window.setSMILES;
  if (typeof original === 'function') window.setSMILES = function () { const result = original.apply(this, arguments); schedule(true); return result; };
  window.addEventListener('load', () => {
    const input = $('smiles-input');
    if (input) input.addEventListener('input', () => schedule(true));
    schedule(true);
    let previous = null;
    setInterval(() => {
      if (document.hidden || typeof getSMILES !== 'function') return;
      const current = getSMILES();
      if (current !== previous) { previous = current; schedule(true); }
    }, 900);
  });
})();
