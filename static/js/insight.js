/* ─── insight.js — Chemistry Guide: know the pocket, know your molecule, predict, dock, learn why ─── */

const GUIDE = { pocket: null, profile: null, match: null, pred: null, run: null, inter: null,
                labels: [], shapes: [], pocketOn: false, linesOn: true, seq: 0, timer: null, currentPose: null };

const GUIDE_CLASSES = {
  acidic:      { color: '#ff4d6d', name: 'Acidic (−)',  tip: 'negatively charged — attracts + groups such as amines' },
  basic:       { color: '#4895ef', name: 'Basic (+)',   tip: 'positively charged — attracts − groups such as acids' },
  polar:       { color: '#2ec4b6', name: 'Polar',       tip: 'H-bond partners for –OH, –NH and C=O groups' },
  hydrophobic: { color: '#f4a261', name: 'Greasy',      tip: 'oily — happy next to carbon-rich parts' },
  aromatic:    { color: '#b388eb', name: 'Aromatic',    tip: 'flat rings — your flat rings can stack on them' },
  cofactor:    { color: '#ffd60a', name: 'Cofactor',    tip: 'a helper molecule that does the chemistry' },
  metal:       { color: '#ffffff', name: 'Metal ion',   tip: 'a metal ion — N/O atoms can bind it' },
  small:       { color: '#8d99ae', name: 'Glycine',     tip: 'tiny and flexible' },
};
const GUIDE_INTER = {
  ionic:       { color: '#ff4d6d', icon: '⚡', name: 'Charge pair',     tip: '+ meets − : the strongest kind of contact' },
  hbond:       { color: '#00e5ff', icon: '💧', name: 'Hydrogen bond',   tip: 'N/O atoms sharing an H — a directional "handshake"' },
  pi:          { color: '#b388eb', icon: '🔵', name: 'Ring stacking',   tip: 'flat rings stacked like coins' },
  hydrophobic: { color: '#f4a261', icon: '🟠', name: 'Greasy contact',  tip: 'oily surfaces pack together and hide from water' },
  metal:       { color: '#ffd60a', icon: '🎯', name: 'Metal binding',   tip: 'a ligand atom bound to a metal ion' },
};
const GUIDE_LEVELS = [[0, 'Lab Intern'], [30, 'Chemist'], [80, 'Medicinal Chemist'], [150, 'Drug Hunter'], [300, 'Pharma Legend']];
const GUIDE_CATS = ['weak', 'medium', 'strong'];
const GUIDE_CAT_LABEL = { weak: '🥱 Weak', medium: '🙂 Medium', strong: '💪 Strong' };

const gh = s => escapeHtml(String(s === undefined || s === null ? '' : s));

// ── XP (kept in this browser only) ───────────────────────────────────────────
function guideXP() { try { return parseInt(localStorage.getItem('medchem_xp') || '0', 10) || 0; } catch (e) { return 0; } }
function guideAddXP(n) { try { localStorage.setItem('medchem_xp', String(guideXP() + n)); } catch (e) {} guideRenderXP(); }
function guideRenderXP() {
  const el = document.getElementById('guide-xp'); if (!el) return;
  const xp = guideXP(); let i = 0;
  GUIDE_LEVELS.forEach((l, k) => { if (xp >= l[0]) i = k; });
  const cur = GUIDE_LEVELS[i], nxt = GUIDE_LEVELS[i + 1];
  const pct = nxt ? Math.round(100 * (xp - cur[0]) / (nxt[0] - cur[0])) : 100;
  el.innerHTML = `🏅 <b>${cur[1]}</b> · ${xp} XP ${nxt ? `<span style="color:var(--text-3)">(${nxt[0] - xp} to ${nxt[1]})</span>` : ''}
    <div style="height:5px;border-radius:3px;background:rgba(255,255,255,.08);margin-top:4px;overflow:hidden">
      <div style="width:${pct}%;height:100%;background:linear-gradient(90deg,#ffd60a,#f72585)"></div></div>`;
}

function guideCategory(e) { return e <= -8.0 ? 'strong' : (e <= -6.5 ? 'medium' : 'weak'); }

// ── loading & rendering ──────────────────────────────────────────────────────
function guideSchedule() { clearTimeout(GUIDE.timer); GUIDE.timer = setTimeout(guideRefresh, 500); }

async function guideRefresh() {
  const body = document.getElementById('guide-sections'); if (!body) return;
  const tid = state.selectedTarget, smiles = (typeof getSMILES === 'function') ? getSMILES() : '';
  if (!tid) { body.innerHTML = '<div style="color:var(--text-3)">Pick a protein target first.</div>'; return; }
  const seq = ++GUIDE.seq;
  const post = (url, data) => apiFetch(url, { method: 'POST', body: JSON.stringify(data) }).catch(e => ({ error: e.message }));
  const [pocket, profile, match] = await Promise.all([
    apiFetch(`/api/insight/pocket/${tid}`).catch(e => ({ error: e.message })),
    smiles ? post('/api/insight/ligand', { smiles }) : Promise.resolve(null),
    smiles ? post('/api/insight/match', { target_id: tid, smiles }) : Promise.resolve(null),
  ]);
  if (seq !== GUIDE.seq) return;
  GUIDE.pocket = pocket; GUIDE.profile = profile; GUIDE.match = match;
  GUIDE.predKey = `${tid}|${smiles}`;
  guideRender();
}

function guideRender() {
  const body = document.getElementById('guide-sections'); if (!body) return;
  body.innerHTML = guidePocketHTML() + guideProfileHTML() + guideMatchHTML() + '<div id="guide-after"></div>';
  if (GUIDE.run && GUIDE.inter) guideRenderAfter();
  guideSyncPocketBtn();
}

function guidePocketHTML() {
  const p = GUIDE.pocket;
  let h = `<div class="gd-sec"><div class="gd-h">① Know the pocket <span>what the protein looks like inside your grid box</span></div>`;
  if (!p || p.error) return h + `<div style="color:var(--text-3)">${gh(p ? p.error : 'Loading…')}</div></div>`;
  const c = p.counts, net = p.net_charge;
  h += `<div><b>${p.n_residues}</b> amino acids line this box. The pocket is <b>${gh(p.character)}</b>.
        Net charge of the lining: <b>${net > 0 ? '+' + net : net}</b>.</div>`;
  h += '<div class="gd-chips" style="margin-top:8px">' + Object.keys(GUIDE_CLASSES).filter(k => c[k]).map(k =>
    `<span class="gd-chip" style="border-color:${GUIDE_CLASSES[k].color};color:${GUIDE_CLASSES[k].color}" title="${gh(GUIDE_CLASSES[k].tip)}">${GUIDE_CLASSES[k].name} × ${c[k]}</span>`).join('') + '</div>';
  h += '<div class="gd-chips" style="margin-top:8px">' + p.residues.slice(0, 32).map(r => {
    const col = (GUIDE_CLASSES[r.cls] || {}).color || '#aaa';
    return `<span class="gd-res" style="background:${col}22;border-color:${col}" title="${gh(r.note || (GUIDE_CLASSES[r.cls] || {}).tip || '')}">${r.note ? '★ ' : ''}${gh(r.res)}${r.resi}</span>`;
  }).join('') + '</div><div class="text-xs" style="color:var(--text-3);margin-top:4px">Closest to the box centre first · ★ = known landmark (hover for why it matters)</div>';
  const notes = p.residues.filter(r => r.note);
  if (notes.length) h += '<ul class="gd-list">' + notes.slice(0, 8).map(r => `<li>★ <b>${gh(r.res)}${r.resi}</b> — ${gh(r.note)}</li>`).join('') + '</ul>';
  if (p.hints.length) h += '<div style="margin-top:8px"><b>What kind of molecule fits?</b></div><ul class="gd-list">' + p.hints.map(x => `<li>${gh(x)}</li>`).join('') + '</ul>';
  h += `<button class="btn btn-sm" id="guide-pocket-btn" onclick="guideTogglePocket()">👁 Show pocket chemistry in 3D</button></div>`;
  return h;
}

function guideProfileHTML() {
  const f = GUIDE.profile;
  let h = `<div class="gd-sec"><div class="gd-h">② Know your molecule <span>the chemical "handles" it can use</span></div>`;
  if (!f) return h + '<div style="color:var(--text-3)">Draw or load a molecule in the designer to see its features.</div></div>';
  if (f.error) return h + `<div style="color:var(--red)">${gh(f.error)}</div></div>`;
  h += '<div class="gd-chips">' + (f.features.length ? f.features.map(x =>
    `<span class="gd-chip" title="${gh(x.tip)}">${x.emoji} ${gh(x.label)}</span>`).join('') : '<span style="color:var(--text-3)">no special groups — a plain hydrocarbon-like molecule</span>') + '</div>';
  h += `<div class="text-xs" style="margin-top:8px;color:var(--text-2)">Size: <b>${f.size_label}</b> (${f.heavy_atoms} heavy atoms) ·
        Greasiness: <b>${f.greasiness}</b> (logP ${f.logp}) · H-bond donors <b>${f.hbd}</b> / acceptors <b>${f.hba}</b> ·
        Flat rings <b>${f.arom_rings}</b> · Floppy bonds <b>${f.rot_bonds}</b></div>${typeof guideRo5HTML === 'function' ? guideRo5HTML(f) : ''}${typeof guideSynthHTML === 'function' ? guideSynthHTML(f) : ''}</div>`;
  return h;
}

function guideMatchHTML() {
  const m = GUIDE.match;
  let h = `<div class="gd-sec"><div class="gd-h">③ Predict, then dock <span>trust your chemistry — then test it</span></div>`;
  if (!m || m.error) return h + `<div style="color:var(--text-3)">${gh(m && m.error ? m.error : 'Once a molecule is loaded, the matchmaker compares it with the pocket.')}</div></div>`;
  const col = m.score >= 68 ? 'var(--green)' : (m.score >= 45 ? '#ffd60a' : 'var(--red)');
  h += `<div><b>Matchmaker score: ${m.score}/100</b> — ${gh(m.verdict)}
        <span class="text-xs" style="color:var(--text-3)">(a rule-of-thumb guess for learning, not a real prediction)</span></div>
        <div style="height:8px;border-radius:4px;background:rgba(255,255,255,.08);margin:6px 0 8px;overflow:hidden">
        <div style="width:${m.score}%;height:100%;background:${col}"></div></div>
        <ul class="gd-list">${m.reasons.map(r => `<li>${r.icon} ${gh(r.text)}</li>`).join('') || '<li>No strong matches or clashes spotted.</li>'}</ul>`;
  const mine = GUIDE.pred && GUIDE.pred.key === GUIDE.predKey ? GUIDE.pred.choice : null;
  h += `<div class="gd-pred"><b>🎯 Your call before docking:</b> how strongly will it bind?
        ${GUIDE_CATS.slice().reverse().map(k => `<button class="btn btn-sm${mine === k ? ' btn-primary' : ''}" onclick="guidePredict('${k}')">${GUIDE_CAT_LABEL[k]}</button>`).join(' ')}
        <div class="text-xs" style="color:var(--text-3);margin-top:4px">${mine ? 'Locked in: <b>' + GUIDE_CAT_LABEL[mine] + '</b> — now click ⚡ Dock Now!' : 'Strong ≤ −8.0 · Medium −6.5 to −8.0 · Weak > −6.5 kcal/mol (more negative = tighter). Correct guess = +10 XP.'}</div></div></div>`;
  return h;
}

function guidePredict(choice) {
  GUIDE.pred = { choice, key: GUIDE.predKey };
  guideRender();
}

// ── after docking ────────────────────────────────────────────────────────────
async function guideAfterDock(r) {
  if (!r || !document.getElementById('guide-sections')) return;
  GUIDE.run = r; GUIDE.inter = null; GUIDE.currentPose = null; GUIDE.verdict = null;
  const actual = guideCategory(r.energy);
  const key = `${r.target_id}|${r.smiles}`;
  if (GUIDE.pred && GUIDE.pred.key === key) {
    const d = Math.abs(GUIDE_CATS.indexOf(GUIDE.pred.choice) - GUIDE_CATS.indexOf(actual));
    const xp = d === 0 ? 10 : (d === 1 ? 4 : 0);
    GUIDE.verdict = { guess: GUIDE.pred.choice, actual, xp };
    if (xp) guideAddXP(xp);
    GUIDE.pred = null;
  } else { GUIDE.verdict = { guess: null, actual, xp: 0 }; }
  await guideLoadInteractions(r.docked_pdb);
}

async function guideLoadInteractions(pdb) {
  const r = GUIDE.run; if (!r || !pdb) return;
  try {
    GUIDE.inter = await apiFetch('/api/insight/interactions', { method: 'POST', body: JSON.stringify({ target_id: r.target_id, docked_pdb: pdb }) });
  } catch (e) { GUIDE.inter = { error: e.message, items: [], counts: {} }; }
  guideRenderAfter();
  if (GUIDE.inter && !GUIDE.inter.error) { await guideWait3D(); guideDrawInteractions(); }
}

function guideRenderAfter() {
  const box = document.getElementById('guide-after'); if (!box || !GUIDE.run) return;
  const r = GUIDE.run, v = GUIDE.verdict, it = GUIDE.inter;
  let h = `<div class="gd-sec"><div class="gd-h">⑤ Why it stuck here <span>what the protein and your molecule are actually doing</span></div>`;
  if (v) {
    h += `<div class="gd-verdict">Result: <b>${r.energy} kcal/mol</b> → <b>${GUIDE_CAT_LABEL[v.actual]}</b> binder. `;
    if (v.guess) h += v.guess === v.actual ? `🎉 You called it (${GUIDE_CAT_LABEL[v.guess]})! <b>+${v.xp} XP</b>`
      : (v.xp ? `You guessed ${GUIDE_CAT_LABEL[v.guess]} — close! <b>+${v.xp} XP</b>` : `You guessed ${GUIDE_CAT_LABEL[v.guess]} — not this time. Look at the contacts below and see what you missed.`);
    else h += '<span style="color:var(--text-3)">Lock in a prediction in step ③ before docking to earn XP.</span>';
    h += '</div>';
  }
  const poses = ((r.engine_metadata || {}).poses || []);
  if (poses.length > 1) h += `<div style="margin:8px 0"><b>Other poses found:</b> ` + poses.map((p, i) =>
    `<button class="btn btn-sm${(GUIDE.poseIdx || 0) === i ? ' btn-primary' : ''}" onclick="guideShowPose(${i})">#${i + 1} · ${p.energy}</button>`).join(' ') +
    `<div class="text-xs" style="color:var(--text-3)">Same molecule, different ways to sit in the pocket. Compare which contacts each one makes.</div></div>`;
  if (!it) h += '<div style="color:var(--text-3)">Analysing contacts…</div>';
  else if (it.error) h += `<div style="color:var(--red)">Could not analyse contacts: ${gh(it.error)}</div>`;
  else {
    const c = it.counts;
    h += `<div style="margin:6px 0">${Object.keys(GUIDE_INTER).filter(k => c[k]).map(k => `${GUIDE_INTER[k].icon} <b>${c[k]}</b> ${GUIDE_INTER[k].name.toLowerCase()}${c[k] > 1 ? 's' : ''}`).join(' · ') || 'No clear contacts found'}
          · <b>${it.buried_pct}%</b> of your molecule is buried in the protein</div>`;
    if (it.items.length) h += '<div style="overflow-x:auto"><table class="gd-table"><tbody>' + it.items.slice(0, 16).map(x => {
      const m = GUIDE_INTER[x.type];
      return `<tr><td title="${gh(m.tip)}">${m.icon} ${m.name}</td><td><b>${gh(x.label)}</b></td><td>${gh(x.detail)}</td></tr>`;
    }).join('') + '</tbody></table></div>';
    h += `<div style="margin-top:8px"><button class="btn btn-sm" onclick="guideToggleLines()">${GUIDE.linesOn ? '🙈 Hide' : '👁 Show'} contact lines in 3D</button></div>
          <div class="text-xs" style="color:var(--text-3);margin-top:6px">Rule of thumb: a snug fit (high % buried) with a few charge pairs / H-bonds <i>and</i> greasy contacts usually scores best.
          Charge pairs and H-bonds make binding <i>specific</i>; greasy contacts add most of the <i>strength</i>. Contacts are estimated from geometry.</div>`;
  }
  box.innerHTML = h + '</div>';
}

async function guideShowPose(i) {
  const r = GUIDE.run, p = ((r && r.engine_metadata) || {}).poses; if (!p || !p[i]) return;
  GUIDE.poseIdx = i; GUIDE.currentPose = p[i].pdb; GUIDE.inter = null; guideRenderAfter();
  await showStructure(r.target_id, { pose: p[i].pdb });
  await guideLoadInteractions(p[i].pdb);
}

// ── 3D helpers ───────────────────────────────────────────────────────────────
function guideClear3D() {
  const v = state.viewer; if (!v) return;
  GUIDE.labels.forEach(l => { try { v.removeLabel(l); } catch (e) {} });
  GUIDE.shapes.forEach(s => { try { v.removeShape(s); } catch (e) {} });
  GUIDE.labels = []; GUIDE.shapes = [];
}

function guideSyncPocketBtn() {
  const b = document.getElementById('guide-pocket-btn');
  if (b) b.textContent = GUIDE.pocketOn ? '🙈 Hide pocket chemistry' : '👁 Show pocket chemistry in 3D';
}

async function guideWait3D() { for (let i = 0; i < 40 && !(state.ligModel && state.proteinModel); i++) await new Promise(r => setTimeout(r, 150)); }

async function guideTogglePocket() {
  if (GUIDE.pocketOn) {
    guideClear3D(); GUIDE.pocketOn = false;
    if (state.proteinModel) { styleProtein('cartoon'); }
    guideSyncPocketBtn(); if (GUIDE.run && GUIDE.inter) guideDrawInteractions(); return;
  }
  const pk = GUIDE.pocket; if (!pk || pk.error) return;
  if (!state.viewer || !state.proteinModel) {
    const run = GUIDE.run && GUIDE.run.target_id === state.selectedTarget ? (GUIDE.currentPose || GUIDE.run.docked_pdb) : null;
    await showStructure(state.selectedTarget, run ? { pose: run } : {});
  }
  const v = state.viewer, p = state.proteinModel; if (!v || !p) return;
  guideClear3D();
  p.setStyle({}, { cartoon: { color: '#7f8fa6', opacity: 0.55, thickness: 0.35 } });
  const groups = {};
  pk.residues.forEach(r => { (groups[`${r.cls}|${r.chain}`] = groups[`${r.cls}|${r.chain}`] || []).push(r.resi); });
  Object.keys(groups).forEach(k => {
    const [cls, chain] = k.split('|'), col = GUIDE_CLASSES[cls].color;
    p.setStyle({ chain, resi: groups[k] }, { cartoon: { color: '#7f8fa6', opacity: 0.55 }, stick: { color: col, radius: 0.2 } });
  });
  pk.residues.filter((r, i) => i < 20 || r.note).forEach(r => {
    GUIDE.labels.push(v.addLabel(`${r.res}${r.resi}`, { position: { x: r.xyz[0], y: r.xyz[1], z: r.xyz[2] }, fontSize: 11, fontColor: 'white',
      backgroundColor: GUIDE_CLASSES[r.cls].color, backgroundOpacity: 0.75, showBackground: true, inFront: true }));
  });
  GUIDE.pocketOn = true; guideSyncPocketBtn();
  const lg = document.getElementById('viewer-legend-guide');
  if (lg) lg.innerHTML = Object.keys(GUIDE_CLASSES).filter(k => pk.counts[k]).map(k => `<span style="color:${GUIDE_CLASSES[k].color}">■</span> ${GUIDE_CLASSES[k].name}`).join(' &nbsp; ');
  if (GUIDE.run && GUIDE.inter) guideDrawInteractions(true);
  v.render();
}

function guideDrawInteractions(keepLabels) {
  const v = state.viewer, p = state.proteinModel, it = GUIDE.inter; if (!v || !it || it.error) return;
  GUIDE.shapes.forEach(s => { try { v.removeShape(s); } catch (e) {} }); GUIDE.shapes = [];
  if (!GUIDE.linesOn) { v.render(); return; }
  it.items.forEach(x => {
    const col = GUIDE_INTER[x.type].color, a = x.lig, b = x.prot;
    const d = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]) || 1, step = 0.45, n = Math.max(1, Math.floor(d / step));
    for (let k = 0; k < n; k += 2) {
      const t0 = k / n, t1 = Math.min(1, (k + 1) / n);
      GUIDE.shapes.push(v.addCylinder({ start: { x: a[0] + (b[0] - a[0]) * t0, y: a[1] + (b[1] - a[1]) * t0, z: a[2] + (b[2] - a[2]) * t0 },
        end: { x: a[0] + (b[0] - a[0]) * t1, y: a[1] + (b[1] - a[1]) * t1, z: a[2] + (b[2] - a[2]) * t1 }, radius: 0.045, color: col, fromCap: 1, toCap: 1 }));
    }
    if (p && !GUIDE.pocketOn) p.setStyle({ chain: x.chain, resi: [x.resi] }, { stick: { color: col, radius: 0.16 } }, true);
  });
  if (!GUIDE.pocketOn) it.items.slice(0, 12).forEach(x => {
    GUIDE.labels.push(v.addLabel(x.label, { position: { x: x.prot[0], y: x.prot[1], z: x.prot[2] }, fontSize: 11, fontColor: 'white',
      backgroundColor: GUIDE_INTER[x.type].color, backgroundOpacity: 0.7, showBackground: true, inFront: true }));
  });
  const lg = document.getElementById('viewer-legend-guide');
  if (lg && !GUIDE.pocketOn) lg.innerHTML = Object.keys(GUIDE_INTER).filter(k => it.counts[k]).map(k => `<span style="color:${GUIDE_INTER[k].color}">┅</span> ${GUIDE_INTER[k].name}`).join(' &nbsp; ');
  v.render();
}

function guideToggleLines() {
  GUIDE.linesOn = !GUIDE.linesOn;
  if (!GUIDE.linesOn) guideClear3D();
  guideDrawInteractions(); guideRenderAfter();
}

// ── hooks into the existing app (runs once, after every other script has loaded) ──
(function guideHooks() {
  const wrap = (name, fn) => { const o = window[name]; if (typeof o === 'function') window[name] = function () { const r = o.apply(this, arguments); try { fn.apply(this, arguments); } catch (e) { console.warn(e); } return r; }; };
  wrap('selectTarget', () => { GUIDE.run = null; GUIDE.inter = null; guideSchedule(); });
  wrap('displayResults', r => guideAfterDock(r));
  wrap('setSMILES', () => guideSchedule());
  wrap('showStructure', () => { GUIDE.pocketOn = false; guideSyncPocketBtn(); });
  const ub = window.updateBox;
  if (typeof ub === 'function') window.updateBox = async function () { const r = await ub.apply(this, arguments); guideSchedule(); return r; };
  window.addEventListener('load', () => {
    const inp = document.getElementById('smiles-input');
    if (inp) inp.addEventListener('input', () => { GUIDE.run = null; guideSchedule(); });
    guideRenderXP(); setTimeout(guideRefresh, 1500); setTimeout(guideRefresh, 5000);
  });
})();
