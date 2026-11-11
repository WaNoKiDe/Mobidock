/* records.js — Cell 49t: the chain of records, the pose tie-break and the record history on the Dock tab */
(function () {
  'use strict';
  if (window.__rec49t) return;
  window.__rec49t = true;

  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const toast = (m, t, i) => { try { showToast(m, t || 'info', i || '🏆'); } catch (e) {} };
  const verdicts = {}, claiming = {}, failedAt = {};
  const keyOf = u => [u.target, u.smiles, u.engine].join('|');
  const parentOf = u => {                                  // the molecule the player started from = the pinned reference molecule
    try { const p = window.MMPin && window.MMPin(); if (p && p.smiles && p.tid === u.target) return p.smiles; } catch (e) {}
    return PL.smiles || u.smiles;
  };

  // The page asks this for the verdict. A clear energy win/loss is final; a tie waits for the pose comparison.
  window.mmOutcome = function (diff, u) {
    const base = diff > PL_TIE ? 'win' : (diff >= -PL_TIE ? 'tie' : 'lose');
    if (base !== 'tie') return base;
    const k = keyOf(u), v = verdicts[k];
    if (!v) return failedAt[k] ? 'tie' : 'pending';
    return (v.decided_by === 'pose' && (v.result === 'win' || v.result === 'lose')) ? v.result : 'tie';
  };

  async function startClaim(u) {
    const k = keyOf(u);
    if (verdicts[k] || claiming[k] || (failedAt[k] && Date.now() - failedAt[k] < 10000)) return;
    claiming[k] = true;
    let v = null;
    try {
      v = await plPost('/api/pipeline/claim', { target_id: u.target, engine: u.engine, parent_smiles: parentOf(u), smiles: u.smiles, player: plPlayer() });
    } catch (e) { failedAt[k] = Date.now(); }
    claiming[k] = false;
    if (v && v.result === 'no-reference') return;           // the reference is still being docked: asked again on the next evaluation
    if (v) {
      verdicts[k] = v;
      if (v.new_record) {
        try { const p = window.MMPin && window.MMPin(); if (p) p.cache = null; } catch (e) {}   // the next refresh picks up the new record
        window.mmRecordsRefresh();
        toast('New record on ' + plProtein() + '! Your molecule (' + u.energy + ' kcal/mol) is now the score to beat for the next player.', 'success', '🏆');
      }
    }
    try { plEvaluate(); } catch (e) {}
    try { plRenderRef(); } catch (e) {}
  }

  // RMSD lines: always shown once the server has judged the docking (not only on ties)
  function rmsdBits(v) {
    const bits = [], rr = v.rmsd_rows || [], best = rr.find(x => x.best);
    if (best) {
      bits.push(best.rmsd != null
        ? '🧬 <b>Pose RMSD</b> vs the best reference (' + esc(best.name) + '): <b>' + esc(best.rmsd) + ' Å</b> — ' + (best.rmsd <= v.same_pose_rmsd ? 'the same binding mode.' : 'a different binding mode.')
        : '🧬 Pose RMSD vs the best reference (' + esc(best.name) + '): not available (too little structure in common, or no pose is stored for that reference).');
    }
    if (v.rmsd_user_native != null) {
      bits.push('📍 Distance to the experimental ligand pose: your pose <b>' + esc(v.rmsd_user_native) + ' Å</b>' +
        (v.rmsd_ref_native != null ? ', best reference pose <b>' + esc(v.rmsd_ref_native) + ' Å</b>.' : '.'));
    }
    if (rr.length > 1) {
      bits.push('Pose RMSD vs each reference: ' + rr.map(x => esc(x.name) + ' (' + esc(x.energy) + ' kcal/mol): ' + (x.rmsd != null ? esc(x.rmsd) + ' Å' : 'n/a')).join(' · '));
    }
    return bits;
  }

  function decorate() {                                    // explain the verdict under the result card
    const u = PL.user; if (!u) return;
    const v = verdicts[keyOf(u)]; if (!v) return;
    const bits = rmsdBits(v);
    if (v.decided_by === 'pose') {
      bits.push('⚖️ The energies tie, so the poses decided: your pose is closer to the experimental ligand pose than the reference pose (or the reverse), by at least 0.5 Å.');
    } else if (v.result === 'tie' && v.rmsd_poses != null) {
      bits.push(v.rmsd_poses <= v.same_pose_rmsd
        ? '⚖️ The energies tie and it is the same binding mode as the reference: a tie, nothing new.'
        : '⚖️ The energies tie and the binding mode differs, but no experimental pose is available to judge which is better: the reference keeps the record.');
    }
    if (v.new_record) bits.push('🏆 <b>New record.</b> Your molecule is stored and is now the score to beat for the next player.');
    if (!bits.length) return;
    ['pl-card', 'pl-result'].forEach(id => {
      const el = document.getElementById(id);
      if (el && !el.querySelector('.mm-claim')) {
        const d = document.createElement('div'); d.className = 'mm-claim text-xs'; d.style.cssText = 'margin-top:8px;color:var(--text-2)';
        d.innerHTML = bits.join('<br>'); el.appendChild(d);
      }
    });
  }

  const origEval = window.plEvaluate;
  window.plEvaluate = function () {
    try {
      const u = PL.user;
      if (u && PL.best && PL.bestEngine === u.engine && u.smiles === PL.smiles && u.target === PL.tid) startClaim(u);
    } catch (e) {}
    const r = origEval.apply(this, arguments);
    try { decorate(); } catch (e) {}
    return r;
  };

  // ── record history under the "Score to beat" panel ──
  const hist = {};
  const histKey = () => {
    const p = window.MMPin && window.MMPin();
    return (p && p.smiles && PL.analysis) ? [PL.tid, plEng(), p.smiles].join('|') : null;
  };
  function paint(el, rows) {
    if (!rows.length) { el.innerHTML = '<div class="text-xs" style="color:var(--text-3)">🏆 No player has beaten this reference yet. The molecule that does becomes the score to beat.</div>'; return; }
    el.innerHTML = '<div class="text-xs"><b>🏆 Record history</b> (newest first)</div>' + rows.slice(0, 5).map(r =>
      '<div class="text-xs" style="color:var(--text-2)">' + esc(r.energy) + ' kcal/mol · ' + esc(r.player) + ' · ' + esc(String(r.created).slice(0, 10)) +
      (r.prev_energy != null ? ' (beat ' + esc(r.prev_energy) + ')' : '') + '</div>').join('');
  }
  function fetchHist(k) {
    const [tid, eng, parent] = [k.split('|')[0], k.split('|')[1], k.split('|').slice(2).join('|')];
    apiFetch('/api/records?target_id=' + encodeURIComponent(tid) + '&engine=' + encodeURIComponent(eng) + '&parent_smiles=' + encodeURIComponent(parent) + '&limit=5')
      .then(d => { hist[k] = d.records || []; const el = document.getElementById('mm-hist'); if (el && histKey() === k) paint(el, hist[k]); })
      .catch(() => { delete hist[k]; });
  }
  window.mmRecordsRefresh = function () { Object.keys(hist).forEach(k => delete hist[k]); try { plRenderRef(); } catch (e) {} };

  const origRef = window.plRenderRef;
  window.plRenderRef = function () {
    origRef.apply(this, arguments);
    try {
      const box = document.getElementById('pl-ref'), k = histKey(); if (!box || !k) return;
      // RMSD column in the reference table: your docked pose vs each reference pose (shown once you have docked)
      const u0 = PL.user, v0 = u0 && verdicts[keyOf(u0)], trs = box.querySelectorAll('table.gd-table tbody tr');
      trs.forEach((tr, i) => {
        const nm = ((tr.querySelector('td:nth-child(2) b') || {}).textContent || '').trim();
        const rr = v0 ? (v0.rmsd_rows || []) : [];
        const m = rr.length === trs.length ? rr[i] : rr.find(x => x.name === nm);
        const td = document.createElement('td'); td.className = 'text-xs'; td.style.cssText = 'white-space:nowrap;color:var(--text-2);text-align:right';
        td.title = 'RMSD between your docked pose and this reference pose';
        td.textContent = !v0 ? '—' : (m && m.rmsd != null ? Number(m.rmsd).toFixed(2) + ' Å' : 'n/a');
        tr.appendChild(td);
      });
      if (trs.length && !v0) {
        const hint = document.createElement('div'); hint.className = 'text-xs'; hint.style.cssText = 'margin-top:4px;color:var(--text-3)';
        hint.textContent = '🧬 The last column (pose RMSD) fills in after you dock your molecule: it compares your pose with each reference pose.';
        const tb = box.querySelector('table.gd-table'); if (tb) tb.insertAdjacentElement('afterend', hint);
      }
      const d = document.createElement('div'); d.id = 'mm-hist'; d.style.cssText = 'margin-top:8px';
      box.appendChild(d);
      if (hist[k] === undefined) { hist[k] = 'loading'; d.innerHTML = '<div class="text-xs" style="color:var(--text-3)">Loading record history…</div>'; fetchHist(k); }
      else if (hist[k] === 'loading') { d.innerHTML = '<div class="text-xs" style="color:var(--text-3)">Loading record history…</div>'; }
      else paint(d, hist[k]);
      const u = PL.user, v = u && verdicts[keyOf(u)];                          // pose RMSD of the player's docking vs each reference
      if (v) {
        const bits = rmsdBits(v);
        if (bits.length) { const r = document.createElement('div'); r.id = 'mm-rmsd'; r.className = 'text-xs'; r.style.cssText = 'margin-top:8px;color:var(--text-2)'; r.innerHTML = bits.join('<br>'); box.appendChild(r); }
      }
    } catch (e) {}
  };
})();
