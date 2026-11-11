/* ─── synth.js — Synthesizability: live meter in the designer, result card after docking, chips elsewhere ─── */

const SYNTH_TIERS = [
  { max: 3.0, key: 'easy',      label: '🟢 Easy to make', colour: '#2af598' },
  { max: 5.0, key: 'moderate',  label: '🟡 Moderate',     colour: '#ffd60a' },
  { max: 6.5, key: 'hard',      label: '🟠 Challenging',  colour: '#f8961e' },
  { max: 99,  key: 'very_hard', label: '🔴 Very hard',    colour: '#ff4d6d' },
];
const SYNTH = { timer: null, seq: 0, last: '' };
const sy = s => escapeHtml(String(s === undefined || s === null ? '' : s));

function synthTierFor(sa) { return SYNTH_TIERS.find(t => sa <= t.max) || SYNTH_TIERS[SYNTH_TIERS.length - 1]; }
function synthMode() { const el = document.getElementById('synth-mode-input'); return el ? el.value : 'standard'; }

// horizontal meter: left = easy (SA 1), right = very hard (SA 10); the marker sits at the molecule's SA score
function synthMeterHTML(sa) {
  const pos = Math.max(0, Math.min(100, ((sa - 1) / 9) * 100));
  return `<div style="position:relative;height:10px;border-radius:5px;margin:8px 0 4px;background:linear-gradient(90deg,#2af598 0%,#ffd60a 33%,#f8961e 55%,#ff4d6d 100%)" title="SA score ${sa} (1 = easy, 10 = very hard)">
      <div style="position:absolute;left:calc(${pos}% - 6px);top:-4px;width:12px;height:18px;border-radius:4px;background:#fff;box-shadow:0 0 0 2px rgba(0,0,0,.55)"></div></div>
    <div class="text-xs" style="display:flex;justify-content:space-between;color:var(--text-3)"><span>1 · easy</span><span>10 · very hard</span></div>`;
}

function synthDriversHTML(s) {
  const e = s.explanation;
  if (!e || e.error) return '<div class="text-xs">Fragment explanation unavailable. Recalculate after restarting the server.</div>';
  const num = x => Number(x).toFixed(4);
  const rows = e.rows.map(r => `<tr><td>${sy(r.id)}</td><td>${r.count}</td><td>${r.matched ? 'Matched' : 'Unseen (default)'}</td><td>${num(r.weight)}</td><td>${num(r.contribution)}</td><td>${r.environments.map(x => `atom ${x.atom}, r=${x.radius}`).join('; ')}</td></tr>`).join('');
  const images = e.rows.filter(r => r.image).map(r => `<figure style="margin:4px;max-width:210px"><img src="${r.image}" alt="Representative Morgan environment ${sy(r.id)}" style="width:200px;background:white;border-radius:8px"><figcaption>ID ${sy(r.id)} · count ${r.count}<br>${r.matched ? 'Matched table entry' : 'Unseen ID'} · contribution ${num(r.contribution)}</figcaption></figure>`).join('');
  return `<div class="text-xs" style="margin-top:8px;color:var(--text-2)"><b>Morgan fragment matches:</b> ${e.matched_unique}/${e.unique} unique IDs found in the SA fragment table (${e.matched_occurrences}/${e.occurrences} occurrences).
  <details style="margin-top:6px"><summary>Why this SA score? View fragments and calculation</summary>
  <p>Overlapping atom-centred environments, radius 0–2 bonds. These are sparse count fingerprints, not similarity matches to a reference drug. “Matched” means the hashed ID exists in RDKit's SA contribution table; it does not prove a synthesis route or identify an exact training molecule.</p>
  <div style="display:flex;flex-wrap:wrap">${images}</div>
  <p>Representative drawings for up to six IDs with the largest absolute weighted contributions. Environments overlap and are not independent chemical building blocks. All IDs are listed below; atom indices are zero-based in the molecule used for scoring.</p>
  <div style="overflow:auto;max-height:320px"><table class="gd-table"><thead><tr><th>Morgan ID</th><th>Count</th><th>Table lookup</th><th>Table weight</th><th>Count × weight / total</th><th>Centres and radii</th></tr></thead><tbody>${rows}</tbody></table></div>
  <p>Fragment term = Σ(count × table weight) / ${e.occurrences} = <b>${num(e.fragment_mean)}</b>. Unseen IDs use weight −4. Higher fragment weights favour easier synthesis; these are raw terms, not final SA points.</p>
  <p>Complexity penalties (raw units): ${Object.entries(e.penalties).map(([k,v])=>`${sy(k)} ${num(v)}`).join(' · ')}.<br>Density/symmetry correction: +${num(e.density_correction)}.</p>
  <p>Raw score = ${num(e.fragment_mean)} − ${num(e.penalty_total)} + ${num(e.density_correction)} = <b>${num(e.raw)}</b>.<br>Mapped SA = 11 − (raw + 5) × 9 / 6.5 = ${num(e.mapped)}. If mapped SA &gt; 8, apply 8 + ln(mapped SA − 8); then clamp to 1–10.<br>Final SA: <b>${s.sa}</b> (lower is easier). Game display: round[100 × (10 − rounded SA) / 9] = <b>${s.score}/100</b>.</p>
  <p>${e.verified ? 'Decomposition verified against the installed RDKit scorer.' : 'Warning: decomposition differs from the installed scorer; the displayed SA remains the scorer output.'} SA is a structural estimate, not a validated synthesis route.</p>
  </details></div>`;
}

function synthUnavailableHTML(s) {
  return `<div class="text-xs" style="color:var(--text-3)">🧪 Synthesizability could not be calculated (${sy(s && s.error || 'scorer unavailable')}). Docking still works; the synthesizability term is left out of the game score.</div>`;
}

// ── live strip under the SMILES box ──────────────────────────────────────────
function synthLiveRender(s, err) {
  const box = document.getElementById('synth-live'); if (!box) return;
  if (err) { box.innerHTML = `<span style="color:var(--text-3)">🧪 Synthesizability: waiting for a valid molecule…</span>`; return; }
  if (!s.available) { box.innerHTML = synthUnavailableHTML(s); return; }
  const gate = synthMode(), limit = gate === 'strict' ? 5.0 : (gate === 'standard' ? 6.5 : null);
  const blocked = limit !== null && s.sa > limit;
  box.innerHTML = `<div style="display:flex;align-items:center;gap:10px;flex-wrap:wrap">
      <b>🧪 Synthesizability</b>
      <span class="gd-chip" style="border-color:${s.colour};color:${s.colour}">${s.tier_label} · SA ${s.sa}</span>
      <span style="color:var(--text-2)"><b>${s.score}</b>/100</span>
      ${blocked ? `<span class="gd-chip" style="border-color:#ff4d6d;color:#ff4d6d">⛔ would be blocked by the '${gate}' gate</span>` : ''}
    </div>${synthMeterHTML(s.sa)}${synthDriversHTML(s)}`;
}

async function synthLive(force) {
  const box = document.getElementById('synth-live'); if (!box) return;
  const smiles = (typeof getSMILES === 'function' ? getSMILES() : '').trim();
  if (!smiles) { box.innerHTML = ''; return; }
  const key = smiles + '|' + synthMode();
  if (!force && key === SYNTH.last) return;
  const seq = ++SYNTH.seq;
  try {
    const s = await apiFetch('/api/synth/analyze', { method: 'POST', body: JSON.stringify({ smiles, mode: synthMode() }) });
    if (seq !== SYNTH.seq) return;
    SYNTH.last = key; synthLiveRender(s, false);
  } catch (e) { if (seq === SYNTH.seq) synthLiveRender(null, true); }
}
function synthSchedule(force) { clearTimeout(SYNTH.timer); SYNTH.timer = setTimeout(() => synthLive(force === true), 450); }

// ── after docking: score cell + full card ─────────────────────────────────────
function synthAfterDock(r) {
  const s = r && r.synth, cell = document.getElementById('score-synth'), card = document.getElementById('synth-card');
  if (cell) {
    if (s && s.available) { cell.textContent = s.score; cell.style.color = s.colour; }
    else cell.textContent = '—';
  }
  const unit = document.getElementById('score-synth-unit');
  if (unit) unit.textContent = s && s.available ? `${s.tier_label.replace(/^\S+\s/, '')} · SA ${s.sa}` : '0 – 100';
  if (!card) return;
  if (!s || !s.available) { card.innerHTML = synthUnavailableHTML(s); return; }
  card.innerHTML = `<div class="card-title" style="margin-top:20px;margin-bottom:0"><span class="dot" style="background:${s.colour}"></span>Synthesizability
      <span class="text-xs" style="margin-left:8px;color:${s.colour}">${s.tier_label} · SA ${s.sa} · gate: ${s.mode === 'off' ? 'off' : (s.passed ? 'PASS' : 'FAIL') + ' · ' + s.mode}</span></div>
    <div class="text-xs" style="color:var(--text-2);margin-top:6px">Easier-to-make molecules earn more game score (up to +2.0). ${s.score}/100 here adds <b>+${(2 * s.norm).toFixed(2)}</b>.</div>
    ${synthMeterHTML(s.sa)}`;
}

// ── chips for the Drug Twin Finder and the Chemistry Guide ───────────────────
function simSynthChip(r) {
  const s = r.synth; if (!s || !s.available) return '';
  const hint = (s.drivers || []).filter(d => d.key !== 'fragments').map(d => d.label).join(', ');
  return `<span class="gd-chip" style="border-color:${s.colour};color:${s.colour}" title="SA score ${s.sa} (1 easy – 10 very hard). Open the Chemistry Guide for Morgan fragment contributions">🧪 ${s.tier_label} · SA ${s.sa}</span>`;
}
function guideSynthHTML(f) {
  const s = f && f.synth; if (!s || !s.available) return '';
  return `<div style="margin-top:8px"><b>🧪 Synthesizability</b> — structural SA estimate
  <div class="gd-chips" style="margin-top:4px"><span class="gd-chip" style="border-color:${s.colour};color:${s.colour}">${s.tier_label} · SA ${s.sa} · ${s.score}/100</span></div>${synthDriversHTML(s)}</div>`;
}

// ── wiring (wraps existing functions instead of editing them, so it survives other patches) ──
(function synthHooks() {
  const wrapFn = (name, after) => {
    const orig = window[name];
    if (typeof orig !== 'function' || orig.__synthWrapped) return;
    const w = function () { const r = orig.apply(this, arguments); try { after.apply(this, arguments); } catch (e) { console.warn('synth hook', name, e); } return r; };
    w.__synthWrapped = true; window[name] = w;
  };
  wrapFn('displayResults', r => synthAfterDock(r));
  wrapFn('setSMILES', () => synthSchedule(true));
  window.addEventListener('load', () => {
    const inp = document.getElementById('smiles-input');
    if (inp) inp.addEventListener('input', () => synthSchedule(false));
    const sel = document.getElementById('synth-mode-input');
    if (sel) sel.addEventListener('change', () => synthSchedule(true));
    synthSchedule(true);
    // the 2-D editor writes the SMILES box without firing an 'input' event, so also watch for silent changes
    let seen = null;
    setInterval(() => {
      if (document.hidden || typeof getSMILES !== 'function') return;
      const cur = getSMILES();
      if (cur !== seen) { if (seen !== null) synthSchedule(false); seen = cur; }
    }, 800);
  });
})();
