
/* ─── pipeline.js — Ligand Pipeline & Challenge: route-specific warnings · synthesizability · closest known drugs · best known docking score to beat ─── */

const PL = { seq: 0, timer: null, analysis: null, smiles: '', tid: '', rows: [], best: null, pending: [], loose: false, same: false,
             rationale: '', bestEngine: '', docking: '', dockSeq: 0, error: '', user: null, lastEnv: '' };
const PL_TIE = 0.3;                                   // kcal/mol: docking noise between two runs is about this big
const pe = s => escapeHtml(String(s === undefined || s === null ? '' : s));
const plCol = v => v >= 60 ? '#2af598' : (v >= 40 ? '#ffd60a' : '#ff4d6d');
const plEng = () => state.selectedEngine || 'vina';
const plExh = () => parseInt((document.getElementById('exhaustiveness-input') || {}).value) || 8;
const plPlayer = () => ((document.getElementById('player-input') || {}).value || '').trim() || 'Anonymous';
const plProtein = () => { const t = state.targets && state.targets[PL.tid]; return t && t.name ? t.name : (PL.tid || 'this protein'); };
const plAuto = () => { const a = document.getElementById('pl-auto'); return !a || a.checked; };
const plPost = (url, body) => apiFetch(url, { method: 'POST', body: JSON.stringify(body) });

// ── analysis (fast, no docking) ──────────────────────────────────────────────
function plSchedule() { clearTimeout(PL.timer); PL.timer = setTimeout(plAnalyze, 900); }

async function plAnalyze() {
  const body = document.getElementById('pl-body'); if (!body) return;
  const smiles = (typeof getSMILES === 'function' ? getSMILES() : '').trim(), tid = state.selectedTarget;
  if (!smiles) { PL.analysis = null; PL.dockSeq++; body.innerHTML = '<div class="gd-sec">Load or draw a ligand and it is screened automatically: route-specific warnings → closest known drugs → the best known docking score for this protein with a molecule like yours.</div>'; return; }
    const seq = ++PL.seq;
  body.innerHTML = '<div class="gd-sec"><span class="spinner" style="width:14px;height:14px;display:inline-block;vertical-align:middle"></span> Screening your ligand…</div>';
  try {
    const a = await plPost('/api/pipeline/analyze', { smiles, target_id: tid || null });
    if (seq !== PL.seq) return;
    PL.dockSeq++;                                     // stop the docking loop that belonged to the previous molecule
    Object.assign(PL, { analysis: a, smiles, tid, rows: [], best: null, pending: [], loose: !!a.plan.loose, same: !!a.plan.same,
                        rationale: a.plan.rationale, bestEngine: '', docking: '', error: '' });
    plRender();
    if (tid) await plRefresh(true);
    else plRenderRef();
  } catch (e) {
    if (seq !== PL.seq) return;
    PL.analysis = null;
    body.innerHTML = `<div class="gd-sec" style="color:var(--text-3)">Waiting for a valid molecule… <span class="text-xs">(${pe(e.message)})</span></div>`;
  }
}

// ── the best known docking score for THIS protein and molecules like THIS ligand ──
async function plRefresh(start) {
  if (!PL.analysis) return;
  const key = PL.smiles + '|' + PL.tid;
  try {
    const d = await plPost('/api/pipeline/best_known', { smiles: PL.smiles, target_id: PL.tid, engine: plEng(), exhaustiveness: plExh(), player: plPlayer() });
    if (key !== PL.smiles + '|' + PL.tid) return;
    Object.assign(PL, { rows: d.rows, best: d.best, pending: d.pending, bestEngine: d.engine, loose: d.loose, same: d.same, rationale: d.rationale, error: '' });
  } catch (e) { PL.error = e.message; }
  plRenderRef(); plRenderBanner(); plEvaluate();
  if (start && plAuto()) plDockPending();
}

// docks the comparison molecules one after another (the server also does one at a time)
async function plDockPending(force) {
  if (PL.docking && !force) return;
  const seq = ++PL.dockSeq, tried = new Set();
  while (seq === PL.dockSeq) {
    const next = PL.pending.find(p => !tried.has(p.smiles));
    if (!next) break;
    tried.add(next.smiles); PL.docking = next.smiles; PL.error = ''; plRenderRef(); plRenderBanner();
    try { await plPost('/api/pipeline/reference_dock', { target_id: PL.tid, smiles: next.smiles, name: next.name, exhaustiveness: plExh(), engine: plEng() }); }
    catch (e) { PL.error = `${next.name}: ${e.message}`; PL.docking = ''; plRenderRef(); break; }
    if (seq !== PL.dockSeq) return;
    PL.docking = ''; await plRefresh(false);
  }
  if (seq === PL.dockSeq) { PL.docking = ''; plRenderRef(); plRenderBanner(); }
}

async function plRedockAll() {
  if (!PL.analysis) return;
  const seq = ++PL.dockSeq;
  for (const c of PL.analysis.plan.candidates) {
    if (seq !== PL.dockSeq) return;
    PL.docking = c.smiles; PL.error = ''; plRenderRef(); plRenderBanner();
    try { await plPost('/api/pipeline/reference_dock', { target_id: PL.tid, smiles: c.smiles, name: c.name, exhaustiveness: plExh(), engine: plEng(), refresh: true }); }
    catch (e) { PL.error = `${c.name}: ${e.message}`; break; }
    await plRefresh(false);
  }
  if (seq === PL.dockSeq) { PL.docking = ''; plRenderRef(); plRenderBanner(); }
}

async function plAddCustom() {
  const smi = (document.getElementById('pl-custom-smiles').value || '').trim(), nm = (document.getElementById('pl-custom-name').value || '').trim();
  if (!smi) { showToast('Paste the SMILES of the molecule to compare against', 'warning', '⚠️'); return; }
  PL.dockSeq++; PL.docking = smi; PL.error = ''; plRenderRef();
  try { await plPost('/api/pipeline/reference_dock', { target_id: PL.tid, smiles: smi, name: nm || 'Added molecule', exhaustiveness: plExh(), engine: plEng() }); }
  catch (e) { PL.error = e.message; }
  PL.docking = ''; await plRefresh(false);
}

// ── rendering ────────────────────────────────────────────────────────────────
const PL_SITE = { site: 'oral' };
function plSetSite(v) { PL_SITE.site = v; plRender(); }
function plDrugHTML(all) {
  const w = all.by_site[PL_SITE.site] || all.by_site.oral;
  const sev = w.blocking.length ? '#ff4d6d' : (w.count ? '#ffd60a' : '#2af598');
  const opts = all.order.map(s => `<option value="${s}"${s === PL_SITE.site ? ' selected' : ''}>${pe(all.labels[s])}${all.summary[s].suitable ? '' : ' ✗'}</option>`).join('');
  const list = (w.groups || []).map(g => `<div style="margin-top:6px"><b>${pe(w.by_group[g][0].group_label)}</b>` +
    w.by_group[g].map(h => `<div class="text-xs" style="margin:3px 0 0 8px;color:${h.severe ? '#ff4d6d' : 'var(--text-2)'}">${h.severe ? '⛔' : '⚠️'} ${pe(h.text)}</div>`).join('') + '</div>').join('');
  return `<div class="gd-sec pl-tile"><div class="gd-h">⚠️ Drug warnings <span>properties that may limit this intended use</span></div>
    <label class="text-xs" style="display:block;margin-bottom:6px">Where must it act?
      <select class="smiles-input" style="width:100%;margin-top:3px" onchange="plSetSite(this.value)">${opts}</select></label>
    <div style="font-weight:700;color:${sev}">${pe(w.verdict)}</div>
    <div class="text-xs" style="color:var(--text-3);margin-top:4px">${w.descriptors_checked} limits checked for this use · ${all.rules} rules in total · ✗ = not suitable for that use</div>
    ${list || '<div class="text-xs" style="color:var(--text-2)">No property crosses its threshold for this use.</div>'}
    <div class="text-xs" style="color:var(--text-3);margin-top:7px">${pe(all.note || 'Rule-of-thumb screen. Delivery, formulation, transporters and clinical data can change the outcome.')}</div></div>`;
}

function plSynthHTML(s) {
  let h = `<div class="gd-sec pl-tile"><div class="gd-h">🧪 Synthesizability <span>can a chemist make it?</span></div>`;
  if (s === null || s === undefined) return h + '<div class="text-xs" style="color:var(--text-3)">The synthesizability check is not installed — run Cell 49k to add it.</div></div>';
  if (!s.available) return h + `<div class="text-xs" style="color:var(--text-3)">Could not be calculated (${pe(s.error || 'scorer unavailable')}).</div></div>`;
  h += `<div style="display:flex;align-items:baseline;gap:8px"><span class="pl-big" style="color:${s.colour}">${s.score}</span><span>/ 100</span><span class="gd-chip" style="border-color:${s.colour};color:${s.colour}">${pe(s.tier_label)} · SA ${s.sa}</span></div>`;
  if (typeof synthMeterHTML === 'function') h += synthMeterHTML(s.sa);
  const d = (s.drivers || []).filter(x => x.key !== 'fragments').slice(0, 2);
  h += `<div class="text-xs" style="color:var(--text-2);margin-top:6px">${d.length ? 'Main difficulty: ' + d.map(x => `${x.icon} ${pe(x.label)}`).join(' · ') : 'No big structural complications.'}</div></div>`;
  return h;
}
function plTwinsHTML(tw) {
  const rows = tw.slice(0, 3).map(t => `<div style="margin:6px 0"><div style="display:flex;justify-content:space-between;gap:8px"><b>${pe(t.name)}</b><span><b>${t.percent}%</b> alike</span></div>
      <div style="height:6px;border-radius:3px;background:rgba(255,255,255,.08);overflow:hidden"><div style="width:${t.percent}%;height:100%;background:linear-gradient(90deg,#4cc9f0,#f72585)"></div></div>
      <div class="text-xs" style="color:var(--text-3)">${pe(t.note)}${t.targets && t.targets.length ? ' · ★ known ligand of ' + pe(t.targets.join(', ')) : ''}</div></div>`).join('');
  const best = tw[0] ? tw[0].percent : 0;
  const msg = best >= 90 ? 'Almost a copy of a known drug.' : (best >= 55 ? 'A relative of known drugs, with differences.' : 'Looks quite different from the known drugs on file.');
  return `<div class="gd-sec pl-tile"><div class="gd-h">🧬 Closest known drugs <span>MACCS · Tanimoto look-alike score</span></div>${rows}
    <div class="text-xs" style="color:var(--text-2)">${msg} The 🧬 Drug Twin Finder (Ligand Browser tab) searches PubChem and ChEMBL for more.</div></div>`;
}



function plRenderRef() {
  const box = document.getElementById('pl-ref'); if (!box) return;
  const a = PL.analysis; if (!a) { box.innerHTML = ''; return; }
  if (!PL.tid) { box.innerHTML = '<div class="gd-sec text-xs">Choose a protein target to compare docking references.</div>'; return; }
  const best = PL.best;
  let h = `<div class="gd-h">🎯 Score to beat <span>best known docking score for ${pe(plProtein())}</span></div>
    <div class="text-xs" style="color:var(--text-2);margin-bottom:6px">${pe(PL.rationale)}</div>`;
  if (best) h += `<div><span class="pl-big" style="color:#4cc9f0">${best.energy}</span> <span>kcal/mol</span>
      <div class="text-xs" style="color:var(--text-3)">${pe(best.name)} · ${best.percent}% alike your ligand · ${pe(plEng())}, exhaustiveness ${plExh()} · runs vary by about ±${PL_TIE}</div></div>`;
  else h += '<div class="text-xs" style="color:var(--text-3)">No docking score yet.</div>';
  h += '<table class="gd-table" style="margin-top:6px"><tbody>' + PL.rows.map(r => {
    const isBest = best && r.smiles === best.smiles;
    const cellv = r.energy !== null && r.energy !== undefined ? `<b>${r.energy}</b>` : (PL.docking === r.smiles ? '<span class="spinner" style="width:11px;height:11px;display:inline-block;vertical-align:middle"></span> docking…' : '<span style="color:var(--text-3)">queued</span>');
    return `<tr><td style="width:16px">${isBest ? '★' : ''}</td><td><b>${pe(r.name)}</b><div class="text-xs" style="color:var(--text-3)">${pe(r.kind_note)}${r.by ? ' · ' + pe(r.by) : ''}</div></td><td>${r.percent}%</td><td style="text-align:right">${cellv}</td></tr>`;
  }).join('') + '</tbody></table>';
  if (PL.error) h += `<div class="text-xs" style="margin-top:6px;color:var(--red)">${pe(PL.error)} <button class="btn btn-sm" onclick="plDockPending(true)">↻ Try again</button></div>`;
  else if (PL.pending.length && !PL.docking) h += `<div style="margin-top:6px"><button class="btn btn-sm btn-primary" onclick="plDockPending(true)">▶ Dock the comparison molecules</button></div>`;
  else if (!PL.docking) h += `<div style="margin-top:6px"><button class="btn btn-sm" onclick="plRedockAll()">↻ Re-dock all</button></div>`;
  h += `<details class="text-xs" style="margin-top:6px;color:var(--text-2)"><summary>Add a molecule to compare against</summary>
      <div style="display:flex;gap:6px;margin-top:6px;flex-wrap:wrap"><input id="pl-custom-name" class="smiles-input" style="width:110px" placeholder="Name">
      <input id="pl-custom-smiles" class="smiles-input" style="flex:1;min-width:140px" placeholder="SMILES"><button class="btn btn-sm" onclick="plAddCustom()">Add and dock</button></div></details>`;
  box.innerHTML = h;
}

function plRenderBanner() {
  const el = document.getElementById('pl-banner'); if (!el) return;
  const a = PL.analysis; if (!a) { el.innerHTML = ''; return; }
  const b = PL.best, left = PL.pending.length;
  if (b) el.innerHTML = `<div class="gd-verdict">🎯 <b>Your mission:</b> dock your ligand and score below <b>${(b.energy - PL_TIE).toFixed(1)} kcal/mol</b> to beat the best known docking score for <b>${pe(plProtein())}</b> with a molecule like yours: <b>${pe(b.name)}</b> (${b.energy}, ${b.percent}% alike). Within ±${PL_TIE} counts as a tie.${left ? ` <span class="text-xs">(${left} more comparison${left > 1 ? 's' : ''} still to dock, so the target may get tougher.)</span>` : ''}</div>`;
  else if (PL.docking || left) el.innerHTML = `<div class="gd-verdict">⏳ <b>Preparing your mission:</b> docking ${left || 1} comparison molecule${left > 1 ? 's' : ''} against ${pe(plProtein())} so you have a score to beat…</div>`;
  else el.innerHTML = '';
}

function plRender() {
  const body = document.getElementById('pl-body'), a = PL.analysis; if (!body || !a) return;
  body.innerHTML = `<div id="pl-banner"></div>
    <div class="pl-grid">${plDrugHTML(a.drug)}${plSynthHTML(a.synth)}${plTwinsHTML(a.twins)}<div class="gd-sec pl-tile" id="pl-ref"></div></div>
    <div id="pl-card"></div>`;
  plRenderRef(); plRenderBanner(); plEvaluate();
}

// ── after the player's own docking: compare with the best known score and grade the ligand ──
function plAfterDock(r) {
  if (!r) return;
  PL.user = { smiles: r.smiles, target: r.target_id, energy: r.energy, engine: r.engine, qed: r.qed };
  if (!PL.analysis || PL.smiles !== r.smiles || PL.tid !== r.target_id) { plAnalyze().then(() => plEvaluate()); return; }
  plRefresh(false);
}


function plAward(n, key) {
  let seen = {}; try { seen = JSON.parse(localStorage.getItem('medchem_pl_awarded') || '{}'); } catch (e) {}
  if (seen[key]) return 0;
  seen[key] = 1;
  try { const ks = Object.keys(seen); if (ks.length > 300) ks.slice(0, 100).forEach(k => delete seen[k]); localStorage.setItem('medchem_pl_awarded', JSON.stringify(seen)); } catch (e) {}
  if (typeof guideAddXP === 'function') guideAddXP(n);
  return n;
}


function plEvaluate() {
  const u = PL.user, a = PL.analysis, best = PL.best;
  const card = document.getElementById('pl-card'), side = document.getElementById('pl-result');
  const put = h => { if (card) card.innerHTML = h; if (side) side.innerHTML = h; };
  if (!u || !a || u.smiles !== PL.smiles || u.target !== PL.tid) return;
  if (!best) { put(`<div class="gd-verdict">⏳ Your score is in (${u.energy} kcal/mol). Waiting for the comparison dockings to finish…</div>`); return; }
  if (PL.bestEngine !== u.engine) {
    put(`<div class="gd-verdict">⚠️ The comparison scores were made with <b>${pe(PL.bestEngine)}</b> but your ligand was docked with <b>${pe(u.engine)}</b>, so they cannot be compared fairly. Select ${pe(u.engine)} in the engine list and the comparison molecules will be docked with it.</div>`);
    return;
  }
  const drug = a.drug.by_site[PL_SITE.site] || a.drug.by_site.oral;
  const diff = +(best.energy - u.energy).toFixed(2);
  const outcome = diff > PL_TIE ? 'win' : (diff >= -PL_TIE ? 'tie' : 'lose');
  const sOk = a.synth && a.synth.available ? a.synth.sa <= 6.5 : null;
  const top = a.twins[0] ? a.twins[0].percent : 0, novel = top < 90;
  let grade;
  if (outcome === 'win' && sOk !== false && novel) grade = 'S';
  else if (outcome === 'win' && sOk !== false) grade = 'A';
  else if (outcome === 'win') grade = 'B';
  else grade = [outcome === 'tie', sOk === true, novel].filter(Boolean).length >= 2 ? 'C' : 'D';
  const key = [u.target, u.smiles, u.engine].join('|');
  const xp = plAward((outcome === 'win' ? 15 : (outcome === 'tie' ? 5 : 0)) + (sOk ? 5 : 0) + (outcome === 'win' && novel ? 5 : 0), key);
  if (xp && outcome === 'win' && typeof showToast === 'function') showToast(`You beat the best known score for ${plProtein()}! +${xp} XP`, 'success', '🏆');
  const title = outcome === 'win' ? `🏆 You beat the best known score (${pe(best.name)}) by ${diff.toFixed(2)} kcal/mol!`
    : (outcome === 'tie' ? `🤝 A tie with ${pe(best.name)} (within ±${PL_TIE} kcal/mol)` : `💪 Not yet — ${pe(best.name)} is ${Math.abs(diff).toFixed(2)} kcal/mol better`);
  const gcol = { S: '#ffd60a', A: '#2af598', B: '#4cc9f0', C: '#f8961e', D: '#ff4d6d' }[grade];
  const row = (ok, label, detail) => `<tr><td style="width:28px">${ok === null ? '➖' : (ok === 'tie' ? '🤝' : (ok ? '✅' : '❌'))}</td><td><b>${label}</b></td><td style="color:var(--text-2)">${detail}</td></tr>`;
  const warningSummary = drug.count ? `${drug.count} threshold warning(s) for ${pe(drug.site_label)}. ${pe(drug.verdict)}` : `No checked thresholds crossed for ${pe(drug.site_label)}.`;
  const tips = outcome === 'lose' ? `<ul class="gd-list"><li>Open the 🧭 Chemistry Guide: which pocket residues are charged or aromatic, and does your molecule have a group to match?</li>
      <li>Check "Why it stuck here" for missing contacts, and flip through the other poses.</li><li>Review the route-specific warnings and synthesizability notes before changing the molecule.</li></ul>` : '';
  put(`<div class="gd-sec" style="border-color:${gcol}">
    <div style="display:flex;align-items:center;gap:14px;flex-wrap:wrap"><div style="font-size:44px;font-weight:900;color:${gcol};line-height:1">${grade}</div>
      <div><div style="font-size:15px;font-weight:700">${title}</div>
      <div class="text-xs" style="color:var(--text-2)">${xp ? `+${xp} XP earned` : 'XP for this ligand was already awarded'} · grade uses docking outcome, synthesizability and novelty; route warnings are shown separately</div></div></div>
    <table class="gd-table" style="margin-top:8px"><tbody>
      ${row(outcome === 'tie' ? 'tie' : (outcome === 'win'), 'Docking vs best known', `yours <b>${u.energy}</b> vs ${pe(best.name)} (${best.percent}% alike) <b>${best.energy}</b> kcal/mol on ${pe(plProtein())} (${diff > 0 ? '−' : '+'}${Math.abs(diff).toFixed(2)})`)}
      ${row(null, 'Drug-use warnings', warningSummary)}
      ${row(sOk, 'Synthesizability', a.synth && a.synth.available ? `${pe(a.synth.tier_label)} · SA ${a.synth.sa}` : 'not available (Cell 49k)')}
      ${row(novel, 'Novelty', `closest known molecule: ${pe(a.twins[0] ? a.twins[0].name : '—')} (${top}% alike)${novel ? '' : ' — nearly a copy'}`)}
    </tbody></table>${tips}</div>`);
}


// ── connection test + source badges for the Drug Twin Finder ─────────────────
function simSourceBadge(r) {
  const chip = (t, c, href) => href ? `<a href="${href}" target="_blank" rel="noopener" class="gd-chip" style="border-color:${c};color:${c};text-decoration:none">${t}</a>` : `<span class="gd-chip" style="border-color:${c};color:${c}">${t}</span>`;
  if (r.source === 'pubchem') return chip(`🌐 PubChem${r.cid ? ' ' + r.cid : ''}`, '#4cc9f0', r.cid ? `https://pubchem.ncbi.nlm.nih.gov/compound/${r.cid}` : '');
  if (r.source === 'chembl') return chip(`🌐 ChEMBL${r.phase === 4 ? ' · approved drug' : (r.phase > 0 ? ' · phase ' + r.phase : '')}`, '#4cc9f0', r.chembl ? `https://www.ebi.ac.uk/chembl/compound_report_card/${r.chembl}/` : '');
  if (r.source === 'known-drug') return chip('💊 known drug', '#2af598');
  if (r.source === 'mine') return chip('✍️ your list', '#8d99ae');
  return r.source ? chip('📚 library', '#8d99ae') : '';
}

async function simDiagnose() {
  const box = document.getElementById('sim-diag'); if (!box) return;
  box.innerHTML = '<span class="spinner" style="width:12px;height:12px;display:inline-block;vertical-align:middle"></span> Testing every route to PubChem and ChEMBL (up to a minute)…';
  try {
    const d = await apiFetch('/api/similarity/diagnose');
    const verdict = !d.any_ok ? '❌ Nothing answered. Restart the local app and retry; otherwise the site may be down or rate-limiting you — wait a minute. Look-alike search still works with the built-in known-drug library.'
      : (d.pubchem_ok && d.chembl_ok ? '✅ Both online databases are reachable.' : (d.pubchem_ok ? '✅ PubChem works. ⚠️ ChEMBL did not answer — PubChem results are still used.' : '⚠️ PubChem did not answer but ChEMBL works — ChEMBL results are used instead.'));
    box.innerHTML = `<div class="gd-sec" style="margin:6px 0"><b>${verdict}</b><table class="gd-table" style="margin-top:6px"><tbody>` +
      d.checks.map(c => `<tr><td style="width:24px">${c.ok ? '✅' : '❌'}</td><td>${pe(c.name)}</td><td style="color:var(--text-2)">${pe(c.detail)}</td><td style="color:var(--text-3)">${c.ms} ms</td></tr>`).join('') + '</tbody></table></div>';
  } catch (e) { box.innerHTML = `<span style="color:var(--red)">The test itself failed: ${pe(e.message)}</span>`; }
}


// ── wiring ───────────────────────────────────────────────────────────────────
(function plHooks() {
  const st = document.createElement('style');
  st.textContent = '.pl-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:10px}.pl-tile{margin-top:0}.pl-big{font-size:34px;font-weight:800;line-height:1}';
  document.head.appendChild(st);
  const wrap = (name, before, after) => { const o = window[name]; if (typeof o !== 'function' || o.__plWrapped) return;
    const w = function () { try { if (before) before.apply(this, arguments); } catch (e) {} const r = o.apply(this, arguments); try { if (after) after.apply(this, arguments); } catch (e) { console.warn('pipeline hook', name, e); } return r; };
    w.__plWrapped = true; window[name] = w; };
  wrap('setSMILES', null, () => plSchedule());
  wrap('selectTarget', null, () => { PL.user = null; plSchedule(); });
  wrap('displayResults', null, r => plAfterDock(r));
  wrap('runDocking', () => { if (PL.analysis && plAuto() && PL.pending.length && !PL.docking) plDockPending(); }, null);
  const ub = window.updateBox;
  if (typeof ub === 'function') window.updateBox = async function () { const r = await ub.apply(this, arguments); plSchedule(); return r; };
  window.addEventListener('load', () => {
    const inp = document.getElementById('smiles-input'); if (inp) inp.addEventListener('input', plSchedule);
    const sp = document.getElementById('score-panel');
    if (sp && !document.getElementById('pl-result')) sp.insertAdjacentHTML('afterbegin', '<div id="pl-result"></div>');
    setTimeout(plAnalyze, 1800);
    setInterval(() => {                                       // engine or exhaustiveness changed → other scores apply
      if (document.hidden || !PL.analysis) return;
      const env = plEng() + '|' + plExh();
      if (PL.lastEnv && env !== PL.lastEnv) { PL.dockSeq++; PL.docking = ''; plRefresh(true); }
      PL.lastEnv = env;
    }, 1000);
  });
})();
