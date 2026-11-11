// ══ CELL 41 · static/js/docking.js ══
/* ─── docking.js — docking — run docking, show scores, engine metadata, property bars ─── */


// ── Docking ───────────────────────────────────────────────────────────────────
async function requestDockResult(engineName, payload) {
  const body = JSON.stringify(payload);
  // All engines run as server-side jobs. DiffDock and GNINA can take minutes,
  // so keeping their work on one long HTTP request makes proxies report
  // "Failed to fetch" even when the docking process is still running.
  const queued = await apiFetch('/api/dock/jobs', { method: 'POST', body });
  if (!queued || !queued.job_id) throw new Error('The server did not return a docking job ID. Restart MobiDock from the updated notebook.');
  const pollUrl = queued.poll_url || `/api/dock/jobs/${encodeURIComponent(queued.job_id)}`;
  let networkRetries = 0;
  const started = Date.now();
  while (true) {
    await new Promise(resolve => setTimeout(resolve, 2000));
    let job;
    try {
      job = await apiFetch(pollUrl);
      networkRetries = 0;
    } catch (error) {
      if (/job not found/i.test(error.message || '')) {
        throw new Error(`The ${engineName} docking job was lost because the server restarted. Please submit it again.`);
      }
      networkRetries++;
      // GNINA can legitimately run for up to 30 minutes. Keep reconnecting
      // through proxy/network interruptions for the full server-side window.
      if (networkRetries > 420) throw new Error(`Lost contact with the ${engineName} docking job after repeated retries: ${error.message}`);
      state.dockingProgress = `Reconnecting to ${engineName} (${networkRetries}/420)…`;
      setLoading(true, state.dockingProgress, 'docking');
      renderDockingRuns();
      await new Promise(resolve => setTimeout(resolve, Math.min(10000, 1500 + networkRetries * 350)));
      continue;
    }
    if (job.status === 'completed' && job.result) return job.result;
    if (job.status === 'failed') throw new Error(job.error || `${engineName} docking job failed.`);
    const seconds = Math.floor((Date.now() - started) / 1000);
    const status = job.status === 'queued' ? 'is queued' : 'is running';
    state.dockingProgress = `${engineName} ${status} (${seconds}s). The job continues on the server while this page stays open.`;
    setLoading(true, state.dockingProgress, 'docking');
    renderDockingRuns();
  }
}

function hideDockingProgress() {
  if (!state.dockingActive) return;
  state.dockingLoaderHidden = true;
  setLoading(false);
  showToast('Progress hidden. Docking continues, and saved poses remain available in the run list.', 'info', 'ℹ️');
}

function pauseDockingAfterCurrent() {
  if (!state.dockingActive) return;
  state.dockingPauseRequested = true;
  const button = document.getElementById('pause-docking');
  if (button) {
    button.disabled = true;
    button.textContent = 'Pause requested…';
  }
  state.dockingProgress = 'Pause requested. The current engine will finish and save its poses before the remaining engines pause.';
  renderDockingRuns();
}

function stopDockingAfterCurrent() {
  if (!state.dockingActive) return;
  state.dockingStopRequested = true;
  const pauseButton = document.getElementById('pause-docking');
  const stopButton = document.getElementById('stop-docking');
  if (pauseButton) pauseButton.disabled = true;
  if (stopButton) {
    stopButton.disabled = true;
    stopButton.textContent = 'Stop requested…';
  }
  state.dockingProgress = 'Stop requested. The current engine will finish and save its poses before the remaining engines stop.';
  renderDockingRuns();
}

function resumeDocking() {
  if (state.dockingActive || !state.dockingResumeBatch) return;
  const resumeBatch = state.dockingResumeBatch;
  runDocking({ resumeBatch });
}

async function runDocking(options = {}) {
  if (state.dockingActive) return;
  const resumeBatch = options.resumeBatch || null;
  const targetId = resumeBatch?.target_id || state.selectedTarget;
  const target = state.targets && state.targets[targetId];
  if (!targetId || !target) {
    showToast('Select a protein target first.', 'error', '⚠️');
    return;
  }

  const ligandSource = document.getElementById('docking-ligand-source')?.value || 'chosen';
  const nativeIndex = document.getElementById('native-ligand-select')?.value;
  const native = resumeBatch ? resumeBatch.native : ligandSource === 'native'
    ? (nativeIndex === '' || nativeIndex == null ? null : state.nativeLigands?.[Number(nativeIndex)])
    : null;
  if (!resumeBatch && ligandSource === 'native' && !native) {
    showToast('Choose an available co-crystal ligand, or select Designer ligand.', 'error', '⚠️');
    return;
  }
  const smiles = resumeBatch?.smiles || native?.smiles || getSMILES();
  if (!smiles) {
    showToast('Load or draw a ligand first.', 'error', '⚠️');
    return;
  }

  const selected = state.engineChoice || state.selectedEngine || 'vina';
  const selectedEngines = resumeBatch
    ? resumeBatch.remaining_engines
    : options.previewFallback
      ? [state.selectedEngine || 'vina']
      : (selected === 'all' ? ['vina', 'gnina', 'diffdock', 'karmadock'] : [selected]);
  const engineNames = selectedEngines.filter(name => {
    const engine = state.engines.find(item => item.name === name);
    return !(name === 'diffdock' && typeof engineRequiresNvidiaGpu === 'function' && engineRequiresNvidiaGpu(engine));
  });
  const requestedPoses = resumeBatch?.requested_pose_count || Math.max(1, Math.min(20, parseInt(document.getElementById('docking-pose-count')?.value, 10) || 10));
  const player = resumeBatch?.player_name || document.getElementById('player-input')?.value.trim() || 'Anonymous';
  const exhaustiveness = resumeBatch?.exhaustiveness || parseInt(document.getElementById('exhaustiveness-input')?.value, 10) || 8;
  const lipinskiMode = resumeBatch?.lipinski_mode || (native ? 'off' : ((document.getElementById('lipinski-mode-input') || {}).value || 'standard'));
  const synthMode = resumeBatch?.synth_mode || (native ? 'off' : ((document.getElementById('synth-mode-input') || {}).value || 'standard'));
  // Native selection preloads the crystallographic box; keep any user edits
  // made to the grid controls when docking that native ligand.
  const box = resumeBatch?.box || readBoxInputs();

  if (native) {
    setSMILES(smiles);
    target.center = box.center;
    target.size = box.size;
  }

  const button = document.getElementById('dock-btn');
  if (button) button.disabled = true;
  state.dockingActive = true;
  state.dockingPauseRequested = false;
  state.dockingStopRequested = false;
  state.dockingPaused = false;
  state.dockingStopped = false;
  state.dockingLoaderHidden = false;
  state.dockingResumeBatch = null;
  state.dockingCurrentEngine = '';
  const pauseButton = document.getElementById('pause-docking');
  const stopButton = document.getElementById('stop-docking');
  if (pauseButton) { pauseButton.disabled = false; pauseButton.textContent = 'Pause after this engine'; }
  if (stopButton) { stopButton.disabled = false; stopButton.textContent = 'Stop after this engine'; }
  const panel = document.getElementById('comparison-results');
  const body = document.getElementById('comparison-table-body');
  state.dockingError = '';
  if (panel && body) {
    panel.hidden = false;
    state.dockingProgress = `Preparing ${native ? (native.ccd_id || 'native ligand') : 'designer ligand'}…`;
    renderDockingRuns();
    panel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  }
  setLoading(true, 'Preparing ligand and docking box…', 'docking');

  const results = [];
  const ligandLabel = resumeBatch?.ligand_label || (native ? `${native.ccd_id || 'Native'} · ${native.name || native.ccd_id || 'PDB ligand'}` : 'Designer ligand');
  const groupKey = resumeBatch?.group_key || `${targetId}\t${(smiles || '').replace(/\s+/g, '')}`;
  const batchContext = {
    target_id: targetId, smiles, native, ligand_label: ligandLabel, group_key: groupKey,
    requested_pose_count: requestedPoses, player_name: player, exhaustiveness,
    lipinski_mode: lipinskiMode, synth_mode: synthMode, box,
  };
  state.dockingActiveContext = batchContext;
  try {
    // Keep the selected/native site box on the server before each run. This also
    // supports older local servers whose /api/dock request model predates box overrides.
    await apiFetch(`/api/target/${encodeURIComponent(targetId)}/box`, {
      method: 'POST', body: JSON.stringify(box),
    });
    target.center = box.center;
    target.size = box.size;

    for (let index = 0; index < engineNames.length; index++) {
      if (state.dockingStopRequested) {
        state.dockingStopRequested = false;
        state.dockingStopped = true;
        state.dockingResumeBatch = null;
        state.dockingProgress = 'Docking stopped. Completed engine poses are saved and downloadable.';
        break;
      }
      if (state.dockingPauseRequested) {
        const remaining = engineNames.slice(index);
        state.dockingPauseRequested = false;
        if (remaining.length) {
          state.dockingResumeBatch = { ...batchContext, remaining_engines: remaining };
          state.dockingPaused = true;
          const nextEngine = state.engines.find(engine => engine.name === remaining[0]);
          state.dockingProgress = `Paused before ${nextEngine?.display_name || remaining[0]}. Completed engine poses are saved and downloadable.`;
        }
        break;
      }
      const engineName = engineNames[index];
      state.dockingCurrentEngine = engineName;
      const engineInfo = state.engines.find(engine => engine.name === engineName) || { name: engineName, display_name: engineName };
      state.dockingProgress = `Running ${engineInfo.display_name || engineName} (${index + 1} of ${engineNames.length})…`;
      renderDockingRuns();
      const diffdockKey = engineName === 'diffdock' && typeof getDiffDockBrowserKey === 'function'
        ? getDiffDockBrowserKey() : '';
      const gpuBlocked = typeof engineRequiresNvidiaGpu === 'function' && engineRequiresNvidiaGpu(engineInfo);
      if (gpuBlocked || (engineInfo.status && engineInfo.status !== 'available' && !diffdockKey)) {
        appendDockingRun({
          target_id: targetId, smiles, ligand_label: ligandLabel, group_key: groupKey,
          engine: engineName, display_name: engineInfo.display_name || engineName,
          status: 'not_configured', message: (engineInfo.problems || []).join('; ') || 'Engine setup is incomplete.',
          pose_count: 0, requested_pose_count: requestedPoses,
        });
        continue;
      }
      setLoading(true, state.dockingProgress, 'docking');
      try {
        const result = await requestDockResult(engineName, {
            target_id: targetId, smiles, player_name: player, exhaustiveness,
            n_poses: requestedPoses, engine: engineName, lipinski_mode: lipinskiMode,
            synth_mode: synthMode, box,
            ...(diffdockKey ? { diffdock_api_key: diffdockKey } : {}),
        });
        results.push(result);
        const metadata = result.engine_metadata || {};
        const poseEnsemble = Array.isArray(metadata.pose_ensemble) && metadata.pose_ensemble.length
          ? metadata.pose_ensemble
          : [{ pdb: result.docked_pdb || result.docked_pdbqt || '', energy: result.energy,
              score_type: metadata.energy_source || 'Binding energy (kcal/mol)' }];
        const runNotes = [];
        const executionModes = [...new Set(poseEnsemble.map(pose => pose.execution_mode || metadata.mode).filter(Boolean))];
        executionModes.forEach(modeValue => {
          const mode = String(modeValue).toLowerCase();
          if (mode.includes('cpu_fallback')) runNotes.push('CPU fallback after GPU request');
          else if (mode.includes('gpu')) runNotes.push('GPU execution');
          else if (mode.includes('cpu')) runNotes.push('CPU execution');
        });
        if (metadata.gpu_fallback_reason) runNotes.push(`GPU fallback: ${metadata.gpu_fallback_reason}`);
        if (metadata.cnn_score != null) runNotes.push(`CNN pose score ${Number(metadata.cnn_score).toFixed(3)}`);
        if (metadata.karma_score != null) runNotes.push(`KarmaDock score ${metadata.karma_score}`);
        if (metadata.box_adjustment?.used_size) {
          runNotes.push(`Vina box adjusted to ${metadata.box_adjustment.used_size.map(value => Number(value).toFixed(1)).join(' × ')} Å per edge to avoid memory failure`);
        }
        const poseRows = poseEnsemble.filter(pose => pose.pdb || pose.text).map((pose, poseIndex) => ({
          target_id: targetId, smiles, ligand_label: ligandLabel, group_key: groupKey,
          engine: result.engine || engineName, display_name: engineInfo.display_name || engineName,
          status: 'completed', energy: pose.energy, confidence: pose.confidence,
          energy_source: pose.score_type || metadata.energy_source || 'Docking score',
          native_score: pose.native_score ?? pose.energy ?? metadata.native_score ?? null,
          native_score_type: pose.native_score_type || pose.score_type || metadata.native_score_type || metadata.energy_source || 'Docking score',
          native_details: pose.native_details || pose.details || metadata.native_details || {},
          execution_mode: pose.execution_mode || metadata.mode || '',
          rmsd: pose.rmsd_to_best_pose ?? pose.rmsd ?? null,
          atoms_compared: pose.atoms_compared ?? null,
          sampling_attempt: pose.attempt ?? pose.sampling_attempt ?? null,
          pose_count: poseEnsemble.length, pose_index: poseIndex + 1,
          requested_pose_count: requestedPoses,
          download: poseIndex === 0 ? (metadata.download || null) : null,
          bundle_url: poseIndex === 0 ? (metadata.download?.url || '') : '',
          pose: pose.pdb || pose.text || '',
          message: poseIndex === 0 ? runNotes.join(' · ') : '',
        }));
        appendDockingRuns(poseRows);
      } catch (error) {
        const notReady = /not set up|not configured|missing|needs setup/i.test(error.message || '');
        appendDockingRun({
          target_id: targetId, smiles, ligand_label: ligandLabel, group_key: groupKey,
          engine: engineName, display_name: engineInfo.display_name || engineName,
          status: notReady ? 'not_configured' : 'failed', message: error.message,
          pose_count: 0, requested_pose_count: requestedPoses,
        });
      }
    }

    // If pause was requested during the last selected engine, there is nothing
    // left to resume; its result has already been saved.
    if (state.dockingPauseRequested) state.dockingPauseRequested = false;
    if (state.dockingStopRequested) {
      state.dockingStopRequested = false;
      state.dockingStopped = true;
      state.dockingResumeBatch = null;
      state.dockingProgress = 'Docking stopped after the current engine. Its completed poses are saved and downloadable.';
    }

    if (state.dockingPaused) {
      showToast('Docking paused. Completed poses are saved and downloadable; resume to run the remaining engines.', 'info', '⏸️');
    } else if (state.dockingStopped) {
      showToast('Docking stopped after the current engine. Completed poses are saved and downloadable.', 'info', '⏹️');
    } else if (results.length) {
      const best = results.reduce((a, b) => Number(b.energy) < Number(a.energy) ? b : a);
      state.lastResult = best;
      displayResults(best);
      showToast(results.length === 1
        ? `Docked with ${best.engine}: ${best.energy} kcal/mol. Results saved.`
        : `${results.length} engine run${results.length === 1 ? '' : 's'} completed; each result is saved separately.`,
      'success', '🎯');
      refreshLeaderboard();
    } else {
      showToast('No selected engine is ready yet. Check each engine status below.', 'warning', '⚠️');
    }
  } catch (error) {
    state.dockingProgress = '';
    state.dockingError = error.message;
    renderDockingRuns();
    showToast(`Docking failed: ${error.message}`, 'error', '❌');
    console.error(error);
  } finally {
    state.dockingActive = false;
    state.dockingCurrentEngine = '';
    state.dockingActiveContext = null;
    if (!state.dockingPaused && !state.dockingStopped) state.dockingProgress = '';
    renderDockingRuns();
    setLoading(false);
    state.dockingLoaderHidden = false;
    if (pauseButton) { pauseButton.disabled = false; pauseButton.textContent = 'Pause after this engine'; }
    if (stopButton) { stopButton.disabled = false; stopButton.textContent = 'Stop after this engine'; }
    if (button) button.disabled = false;
  }
}

// ── Display results ───────────────────────────────────────────────────────────
function displayResults(r) {
  const panel = document.getElementById('score-panel');
  panel.classList.add('visible');
  panel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });

  animateValue(document.getElementById('score-energy'), 0, r.energy, 2);

  // Engine badge
  const engineBadge = document.getElementById('engine-badge');
  const engInfo = state.engines.find(e => e.name === r.engine) || {};
  engineBadge.innerHTML = `${engInfo.icon || '⚗️'} ${engInfo.display_name || r.engine}`;
  engineBadge.style.borderColor = engInfo.color || 'var(--accent)';
  engineBadge.style.color = engInfo.color || 'var(--accent)';
  engineBadge.style.background = `${engInfo.color || 'var(--accent)'}22`;

  const mockBadge = document.getElementById('mock-badge');
  if (mockBadge) mockBadge.style.display = 'none';   // results are always real now

  renderProps(r.props, r.qed);
  /* BADGE_MARK */
  const lb = document.getElementById('lipinski1-badge');
  if (lb) {
    const l = r.props && r.props.lipinski1;
    lb.textContent = !l ? '' : (l.mode === 'off' ? 'Lipinski gate: off' : `Lipinski gate: ${l.passed ? 'PASS' : 'FAIL'} · ${l.mode} · ${l.violations} miss${l.violations === 1 ? '' : 'es'}`);
    lb.style.color = !l || l.mode === 'off' ? 'var(--text-3)' : (l.passed ? 'var(--green)' : 'var(--red)');
  }
  renderEngineMeta(r);
  // the server returns the pose as PDB text for every engine
  showStructure(r.target_id, { pose: r.docked_pdb || r.docked_pdbqt || '' });
  loadLiterature(r.target_id);
}

function renderEngineMeta(r) {
  const el = document.getElementById('engine-meta');
  if (!el) return;
  const m = r.engine_metadata || {};
  const bits = [];
  const f = (v, d = 2) => (v === null || v === undefined || isNaN(v)) ? null : Number(v).toFixed(d);
  if (m.confidence !== undefined && m.confidence !== null) { /* shown below */ }
  if (m.seed !== undefined && m.seed !== null) bits.push(`random seed: <b>${m.seed}</b>`);
  if (m.mode) bits.push(`mode: <b>${escapeHtml(m.mode)}</b>`);
  if (f(m.cnn_score, 3)) bits.push(`CNN pose score: <b>${f(m.cnn_score, 3)}</b>`);
  if (f(m.cnn_affinity)) bits.push(`CNN affinity: <b>${f(m.cnn_affinity)}</b> pK`);
  if (f(m.confidence, 2)) bits.push(`DiffDock confidence: <b>${f(m.confidence, 2)}</b>`);
  if (f(m.karma_score)) bits.push(`KarmaDock score: <b>${f(m.karma_score)}</b>`);
  if (m.energy_source) bits.push('energy = ' + escapeHtml(m.energy_source));
  if (m.pose_inside_grid_box === false) bits.push('<span style="color:var(--accent-4)">⚠ pose centre lies outside the grid box</span>');
  el.innerHTML = bits.join(' &nbsp;·&nbsp; ');
  el.style.display = bits.length ? 'block' : 'none';
}

function appendDockingRun(row) { appendDockingRuns([row]); }

function dockingEnergy(row) {
  if (row.energy === null || row.energy === undefined || String(row.energy).trim() === '') return null;
  const value = Number(row.energy);
  return Number.isFinite(value) ? value : null;
}

function nativeScoreLabel(key) {
  const labels = {
    minimizedAffinity: 'GNINA minimized affinity (kcal/mol)',
    CNNscore: 'GNINA CNN pose score (unitless)',
    CNNaffinity: 'GNINA CNN affinity (pKd)',
    CNN_VS: 'GNINA CNN virtual screening score',
    cnn_score: 'GNINA CNN pose score (unitless)',
    cnn_affinity: 'GNINA CNN affinity (pKd)',
    karma_score: 'KarmaDock model score (native units)',
    karma_score_ff_corrected: 'KarmaDock force-field-corrected score (native units)',
    karma_score_ff: 'KarmaDock force-field-corrected score (native units)',
    karma_score_aligned: 'KarmaDock aligned score (native units)',
  };
  return labels[key] || String(key).replace(/[_-]+/g, ' ').replace(/\b\w/g, letter => letter.toUpperCase());
}

function nativeScoreDetailsHtml(details) {
  let values = details || {};
  if (typeof values === 'string') {
    try { values = JSON.parse(values); } catch (_) { return escapeHtml(values); }
  }
  if (!values || typeof values !== 'object' || Array.isArray(values)) return '—';
  const entries = Object.entries(values).filter(([, value]) => value !== null && value !== undefined && value !== '');
  if (!entries.length) return '<span class="submetric">No additional native score components</span>';
  return `<div class="docking-native-scores">${entries.map(([key, value]) => {
    const numeric = typeof value === 'number' ? value : (String(value).trim() !== '' ? Number(value) : NaN);
    const shown = Number.isFinite(numeric) ? String(Number(numeric.toFixed(4))) : String(value);
    return `<div><span>${escapeHtml(nativeScoreLabel(key))}</span><b>${escapeHtml(shown)}</b></div>`;
  }).join('')}</div>`;
}

function dockingExecutionLabel(mode, engine) {
  const value = String(mode || '').toLowerCase();
  if (value.includes('cpu_fallback')) return 'CPU (GPU fallback)';
  if (value.includes('gpu')) return 'GPU';
  if (value.includes('cpu') || engine === 'vina') return 'CPU';
  if (engine === 'diffdock') return 'Remote DiffDock';
  return value || 'Not reported';
}

function appendDockingRuns(rows) {
  state.dockingRuns = state.dockingRuns || [];
  state.dockingRmsdBest = state.dockingRmsdBest || {};
  state.dockingRmsdEpoch = state.dockingRmsdEpoch || {};
  state.dockingRmsdPending = state.dockingRmsdPending || [];
  const batchId = ++state.dockingBatchSequence;
  const groups = new Set();
  rows.forEach(row => {
    row.run_id = ++state.dockingRunSequence;
    row.batch_id = batchId;
    row.rmsd = null;
    row.atoms_compared = null;
    state.dockingRuns.push(row);
    if (row.status === 'completed' && row.pose) groups.add(row.group_key);
  });
  renderDockingRuns();
  groups.forEach(groupKey => {
    const pending = updateDockingRmsds(groupKey);
    state.dockingRmsdPending.push(pending);
    const removePending = () => {
      state.dockingRmsdPending = (state.dockingRmsdPending || []).filter(item => item !== pending);
    };
    pending.then(removePending, removePending);
  });
}

async function updateDockingRmsds(groupKey) {
  const rows = (state.dockingRuns || []).filter(row => row.group_key === groupKey && row.status === 'completed' && row.pose);
  if (!rows.length) return;
  const scored = rows.filter(row => dockingEnergy(row) !== null);
  if (!scored.length) return;
  const best = scored.reduce((a, b) => dockingEnergy(b) < dockingEnergy(a) ? b : a);
  const bestChanged = state.dockingRmsdBest[groupKey] !== best.run_id;
  if (bestChanged) {
    state.dockingRmsdBest[groupKey] = best.run_id;
    state.dockingRmsdEpoch[groupKey] = (state.dockingRmsdEpoch[groupKey] || 0) + 1;
    rows.filter(row => row !== best).forEach(row => { row.rmsd = null; row.atoms_compared = null; });
  }
  const epoch = state.dockingRmsdEpoch[groupKey] || 0;
  best.rmsd = 0;
  best.atoms_compared = null;
  await Promise.all(rows.filter(row => row !== best && row.rmsd === null).map(async row => {
    try {
      const result = await apiFetch('/api/dock/rmsd', {
        method: 'POST', body: JSON.stringify({ pose_a: best.pose, pose_b: row.pose }),
      });
      if (state.dockingRmsdEpoch[groupKey] !== epoch) return;
      row.rmsd = result.available ? result.rmsd_angstrom : 'n/a';
      row.atoms_compared = result.atoms_compared || null;
    } catch (error) {
      if (state.dockingRmsdEpoch[groupKey] !== epoch) return;
      row.rmsd = 'n/a';
      row.rmsd_message = error.message;
    }
  }));
  if (state.dockingRmsdEpoch[groupKey] === epoch) renderDockingRuns();
}

function toggleDockingBatch(batchId) {
  state.dockingExpandedBatches = state.dockingExpandedBatches || new Set();
  if (state.dockingExpandedBatches.has(Number(batchId))) state.dockingExpandedBatches.delete(Number(batchId));
  else state.dockingExpandedBatches.add(Number(batchId));
  renderDockingRuns();
}

function renderDockingRuns() {
  const panel = document.getElementById('comparison-results');
  const body = document.getElementById('comparison-table-body');
  const note = document.getElementById('comparison-note');
  if (!panel || !body) return;
  const esc = value => escapeHtml(String(value ?? ''));
  panel.hidden = false;
  const rows = state.dockingRuns || [];
  const batches = new Map();
  rows.forEach(row => {
    const key = Number(row.batch_id || row.run_id);
    if (!batches.has(key)) batches.set(key, []);
    batches.get(key).push(row);
  });
  const latest = rows[rows.length - 1];
  const latestGroupRows = latest
    ? rows.filter(row => row.group_key === latest.group_key && row.status === 'completed' && dockingEnergy(row) !== null)
    : [];
  const bestReference = latestGroupRows.reduce((best, row) => !best || dockingEnergy(row) < dockingEnergy(best) ? row : best, null);
  const downloadButton = document.getElementById('download-docking-table');
  if (downloadButton) downloadButton.disabled = rows.length === 0 && !state.dockingActive;
  const resumeButton = document.getElementById('resume-docking');
  if (resumeButton) {
    resumeButton.hidden = !state.dockingPaused;
    resumeButton.disabled = state.dockingActive;
  }
  const compatibility = state.serverCompatibilityWarning ? `${state.serverCompatibilityWarning} · ` : '';
  const poseRowsTotal = rows.filter(row => row.status === 'completed' && row.pose).length;
  note.textContent = compatibility + (rows.length
    ? `${batches.size} docking runs · ${poseRowsTotal} pose${poseRowsTotal === 1 ? '' : 's'} saved · newest: ${latest?.ligand_label || 'Ligand'} on ${latest?.target_id || 'target'}` +
      (bestReference ? ` · RMSD reference: ${bestReference.engine} pose ${bestReference.pose_index || 1} at ${dockingEnergy(bestReference).toFixed(2)} kcal/mol` : ' · waiting for a binding-energy pose')
    : 'Select one engine or All four engines, then Dock Now. Completed results stay in this list.');
  state.dockingExpandedBatches = state.dockingExpandedBatches || new Set();
  const rmsdLabel = row => row.rmsd === 0 ? '0.00 Å (best)'
    : typeof row.rmsd === 'number' ? `${row.rmsd.toFixed(2)} Å`
    : row.rmsd === 'n/a' ? 'n/a'
    : !rows.some(item => item.group_key === row.group_key && dockingEnergy(item) !== null) ? 'n/a · no energy reference'
    : row.status === 'completed' ? 'Calculating…' : '—';
  const html = [...batches.entries()].map(([batchId, runRows]) => {
    const poses = runRows.filter(row => row.status === 'completed' && row.pose);
    const energyRows = poses.filter(row => dockingEnergy(row) !== null);
    const bestPose = energyRows.length
      ? energyRows.reduce((a, b) => dockingEnergy(b) < dockingEnergy(a) ? b : a)
      : poses.reduce((a, b) => Number(b.confidence || -Infinity) > Number(a.confidence || -Infinity) ? b : a, poses[0]);
    const row = bestPose || runRows[0];
    const completed = poses.length > 0;
    const score = !bestPose ? '—' : bestPose.energy == null
      ? (bestPose.confidence == null ? '—' : `confidence ${Number(bestPose.confidence).toFixed(3)}`)
      : `${Number(bestPose.energy).toFixed(2)} kcal/mol`;
    const requested = Number(row.requested_pose_count || 0);
    const poseCount = poses.length;
    const returned = completed ? `${poseCount} pose${poseCount === 1 ? '' : 's'}${requested && requested !== poseCount ? ` · requested ${requested}` : ''}` : '—';
    const statusText = completed ? 'Completed · saved'
      : row.status === 'not_configured' ? 'Setup needed' : 'Failed';
    const message = row.message ? `<span class=\"submetric\">${esc(row.message)}</span>` : '';
    const rmsd = bestPose ? rmsdLabel(bestPose) : '—';
    const zipRow = runRows.find(item => item.download?.url);
    const zip = zipRow ? `<a class=\"btn btn-sm\" href=\"${esc(zipRow.download.url)}\" onclick=\"event.stopPropagation()\">Download run</a>` : (completed ? 'Saved' : '—');
    const rowAttrs = completed
      ? `class=\"docking-run-summary-row\" role=\"button\" tabindex=\"0\" onclick=\"toggleDockingBatch(${batchId})\" onkeydown=\"if(event.target===event.currentTarget&&(event.key==='Enter'||event.key===' ')){event.preventDefault();toggleDockingBatch(${batchId})}\"`
      : '';
    const toggle = completed
      ? `<button class=\"btn btn-sm docking-run-toggle\" type=\"button\" onclick=\"event.stopPropagation();toggleDockingBatch(${batchId})\">${state.dockingExpandedBatches.has(batchId) ? 'Hide poses' : `Show ${poseCount} poses`}</button>`
      : '';
    const summary = `<tr ${rowAttrs}><td><b>${esc(row.display_name || row.engine)}</b><span class=\"submetric\">${esc(row.engine)} · run ${batchId}</span><span class=\"submetric\">${esc(row.target_id)} · ${esc(row.ligand_label)}</span></td><td><span class=\"comparison-status ${esc(row.status)}\">${statusText}</span>${message}</td><td>${score}<span class=\"submetric\">${esc(bestPose?.energy_source || row.energy_source || '')}</span></td><td title=\"RMSD in the protein coordinate frame against the pose with the lowest binding energy for this ligand and target\">${rmsd}${bestPose?.atoms_compared ? `<span class=\"submetric\">${bestPose.atoms_compared} atoms</span>` : ''}</td><td>${returned}</td><td>${zip}</td><td>${toggle}</td></tr>`;
    if (!completed || !state.dockingExpandedBatches.has(batchId)) return summary;
    const detailRows = poses.map(pose => {
      const poseScore = pose.energy == null
        ? (pose.confidence == null ? '—' : `confidence ${Number(pose.confidence).toFixed(3)}`)
        : `${Number(pose.energy).toFixed(2)} kcal/mol`;
      const rmsdDetail = pose.atoms_compared ? `<span class=\"submetric\">${pose.atoms_compared} atoms</span>` : pose.rmsd_message ? `<span class=\"submetric\">${esc(pose.rmsd_message)}</span>` : '';
      const nativeScore = pose.native_score == null ? '—' : (Number.isFinite(Number(pose.native_score)) ? Number(pose.native_score).toFixed(4) : esc(pose.native_score));
      const nativeType = pose.native_score_type ? `<span class=\"submetric\">${esc(pose.native_score_type)}</span>` : '';
      return `<tr><td>Pose ${pose.pose_index || '—'}</td><td>${poseScore}<span class=\"submetric\">${esc(pose.energy_source || '')}</span></td><td><b>${nativeScore}</b>${nativeType}</td><td>${nativeScoreDetailsHtml(pose.native_details)}</td><td>${esc(dockingExecutionLabel(pose.execution_mode, pose.engine))}</td><td>${rmsdLabel(pose)}${rmsdDetail}</td><td>${pose.confidence == null ? '—' : Number(pose.confidence).toFixed(3)}</td><td><button class=\"btn btn-sm comparison-view-button\" type=\"button\" onclick=\"viewDockingRun(${pose.run_id})\">View pose</button></td></tr>`;
    }).join('');
    return summary + `<tr class=\"docking-run-details-row\"><td colspan=\"7\"><div class=\"docking-run-details\"><table class=\"docking-pose-table\"><thead><tr><th>Pose</th><th>Comparable binding energy / score</th><th>Native engine score</th><th>Native score details</th><th>Execution device</th><th>RMSD to best-energy pose</th><th>Confidence</th><th></th></tr></thead><tbody>${detailRows}</tbody></table></div></td></tr>`;
  }).join('');
  const progressAction = state.dockingActive && state.loadingKind !== 'general'
    ? `${state.dockingPauseRequested || state.dockingStopRequested ? '' : `<button class=\"btn btn-sm docking-pause-inline\" type=\"button\" onclick=\"pauseDockingAfterCurrent()\">Pause after this engine</button><button class=\"btn btn-sm docking-stop-inline\" type=\"button\" onclick=\"stopDockingAfterCurrent()\">Stop after this engine</button>`}${state.dockingPauseRequested ? '<span>Pause requested…</span>' : ''}${state.dockingStopRequested ? '<span>Stop requested…</span>' : ''}`
    : '';
  const progress = state.dockingProgress ? `<tr class=\"docking-progress-row\"><td colspan=\"7\"><div class=\"docking-progress-content\"><span>${esc(state.dockingProgress)}</span>${progressAction}</div></td></tr>` : '';
  const error = state.dockingError ? `<tr><td colspan=\"7\" style=\"color:var(--red)\">Docking could not start: ${esc(state.dockingError)}</td></tr>` : '';
  body.innerHTML = html + progress + error || `<tr><td colspan=\"7\">No docking runs yet.</td></tr>`;
}
async function downloadDockingTable() {
  if (!(state.dockingRuns || []).length && !state.dockingActive) return;
  const button = document.getElementById('download-docking-table');
  if (button) { button.disabled = true; button.textContent = 'Preparing Excel…'; }
  try {
    while ((state.dockingRmsdPending || []).length) {
      await Promise.allSettled([...(state.dockingRmsdPending || [])]);
    }
    const rows = [...(state.dockingRuns || [])];
    if (state.dockingActive) {
      const context = state.dockingActiveContext || {};
      const currentEngine = state.dockingCurrentEngine || '';
      const engineInfo = state.engines.find(engine => engine.name === currentEngine) || {};
      rows.push({
        target_id: context.target_id || state.selectedTarget || '',
        ligand_label: context.ligand_label || 'Docking in progress', smiles: context.smiles || '',
        display_name: engineInfo.display_name || (currentEngine ? currentEngine : 'Preparing docking'),
        engine: currentEngine, batch_id: state.dockingBatchSequence + 1,
        status: 'running', pose_index: '', pose_count: 0,
        requested_pose_count: context.requested_pose_count || '',
        energy: '', energy_source: '', native_score: '', native_score_type: '', native_details: '',
        confidence: '', rmsd: '', atoms_compared: '', sampling_attempt: '', bundle_url: '',
        message: state.dockingProgress || 'Docking is in progress; no pose has returned yet.',
      });
    }
    const columns = [
      'target_id', 'ligand_label', 'smiles', 'display_name', 'engine', 'batch_id', 'status',
      'pose_index', 'pose_count', 'requested_pose_count', 'energy', 'energy_source',
      'native_score', 'native_score_type', 'native_details', 'confidence', 'rmsd',
      'atoms_compared', 'sampling_attempt', 'bundle_url', 'message',
    ];
    const payloadRows = rows.map(row => {
      const result = {};
      columns.forEach(key => {
        const value = row[key];
        result[key] = value && typeof value === 'object' ? JSON.stringify(value) : value;
      });
      return result;
    });
    const response = await fetch('/api/docking-results/export.xlsx', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ rows: payloadRows }),
    });
    if (!response.ok) {
      let message = `Excel export failed (${response.status})`;
      try { const body = await response.json(); message = body.detail || body.error || message; } catch (_) {}
      throw new Error(message);
    }
    const blob = await response.blob();
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    const stamp = new Date().toISOString().replace(/[:.]/g, '-');
    link.href = url;
    link.download = `mobidock-docking-results-${stamp}.xlsx`;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  } catch (error) {
    showToast(error.message || 'Could not download the Excel workbook.', 'error', '❌');
  } finally {
    if (button) { button.textContent = 'Download Excel (.xlsx)'; button.disabled = !(state.dockingRuns || []).length && !state.dockingActive; }
  }
}

function viewDockingRun(runId) {
  const row = (state.dockingRuns || []).find(item => item.run_id === runId);
  if (!row || !row.pose) return;
  if (state.selectedTarget !== row.target_id && typeof selectTarget === 'function') selectTarget(row.target_id);
  else state.selectedTarget = row.target_id;
  showStructure(row.target_id, { pose: row.pose });
  const score = row.energy == null
    ? (row.confidence == null ? 'pose' : `confidence ${Number(row.confidence).toFixed(3)}`)
    : `${Number(row.energy).toFixed(2)} kcal/mol`;
  showToast(`${row.display_name || row.engine}: ${score}.`, 'success', '🧬');
}

function renderProps(props, qed) {
  const defs = [
    { key: 'mw',        label: 'Mol. Weight', unit: 'Da',  max: 500,  good: v => v <= 500 },
    { key: 'logp',      label: 'LogP',         unit: '',   max: 5,    good: v => v <= 5 },
    { key: 'hbd',       label: 'H-Bond Donors',unit: '',   max: 5,    good: v => v <= 5 },
    { key: 'hba',       label: 'H-Bond Acc.',  unit: '',   max: 10,   good: v => v <= 10 },
    { key: 'tpsa',      label: 'TPSA',         unit: 'Å²', max: 140,  good: v => v <= 140 },
    { key: 'rot_bonds', label: 'Rot. Bonds',   unit: '',   max: 10,   good: v => v <= 10 },
  ];
  const grid = document.getElementById('props-grid');
  if (!grid) return;
  grid.innerHTML = '';
  defs.forEach(d => {
    const val = props[d.key] ?? 0;
    const pct = Math.min(100, (val / d.max) * 100);
    const ok = d.good(val);
    const color = ok ? 'var(--green)' : 'var(--red)';
    const item = document.createElement('div');
    item.className = 'prop-item';
    item.innerHTML = `
      <div class="prop-label">${d.label}</div>
      <div class="prop-value" style="color:${color}">${val}${d.unit ? ' <span style="font-size:11px;color:var(--text-3)">' + d.unit + '</span>' : ''}</div>
      <div class="prop-bar"><div class="prop-bar-fill" style="width:${pct}%;background:${ok ? 'linear-gradient(90deg,var(--green),#2af598)' : 'linear-gradient(90deg,var(--red),#f4895f)'}"></div></div>
    `;
    grid.appendChild(item);
  });
}

// === LIGAND_PDB_PICKER_V2 ===
/* Search RCSB for the loaded ligand and allow importing a matching receptor.
   This intentionally has no reference score, score comparison, or score-all flow. */
(function () {
  'use strict';
  if (window.__ligandPdbPickerV2) return;
  window.__ligandPdbPickerV2 = true;
  let timer = null, sequence = 0, lastQuery = '', activeSearch = null;
  const statePdb = { query: '', hits: [], loading: false, error: '' };
  const esc = value => escapeHtml(String(value == null ? '' : value));

  function panel() {
    let box = document.getElementById('ligand-pdb-results');
    if (box) return box;
    const grid = document.getElementById('target-grid');
    if (!grid || !grid.parentElement) return null;
    box = document.createElement('div');
    box.id = 'ligand-pdb-results';
    box.style.cssText = 'margin:10px 0 14px;padding:12px;border:1px solid var(--border);border-radius:10px;background:var(--bg-1)';
    grid.parentElement.insertBefore(box, grid);
    return box;
  }

  function render() {
    const box = panel();
    if (!box) return;
    const searchButton = statePdb.query
      ? `<button class="btn btn-sm" type="button" data-ligand-pdb-search>${statePdb.loading ? 'Cancel search' : 'Find PDBs with this native ligand'}</button>`
      : '';
    if (statePdb.loading) {
      box.innerHTML = `<div class="text-xs" style="color:var(--text-2);display:flex;align-items:center;justify-content:space-between;gap:10px;flex-wrap:wrap"><span>Searching RCSB for structures containing this ligand…</span>${searchButton}</div>`;
      box.querySelector('[data-ligand-pdb-search]')?.addEventListener('click', cancelSearch);
      return;
    }
    if (!statePdb.query) {
      box.innerHTML = '<div class="text-xs" style="color:var(--text-3)">Load a ligand to find matching RCSB PDB structures. Choose a result to load its protein as the docking target.</div>';
      return;
    }
    if (!statePdb.hits.length) {
      box.innerHTML = `<div class="text-xs" style="color:var(--text-2);display:flex;align-items:center;justify-content:space-between;gap:10px;flex-wrap:wrap"><span><b>RCSB structures for this ligand</b><br>${statePdb.error ? esc(statePdb.error) : 'No exact co-crystal match was found.'}</span>${searchButton}</div>`;
      box.querySelector('[data-ligand-pdb-search]')?.addEventListener('click', () => search(statePdb.query));
      return;
    }
    box.innerHTML = `<div class="text-xs" style="margin-bottom:8px"><b>RCSB structures containing the loaded ligand</b><br><span style="color:var(--text-3)">Choose a PDB ID to load its protein. The native ligand is removed and the docking box is centred on its binding site.</span> ${searchButton}</div>` +
      statePdb.hits.map((hit, i) => `<div style="display:flex;gap:8px;justify-content:space-between;align-items:center;flex-wrap:wrap;padding:8px;margin-top:6px;border:1px solid var(--border);border-radius:8px"><div class="text-xs"><b>${esc(hit.pdb_id)}</b> · ligand ${esc(hit.ccd_id)}${hit.resolution != null ? ` · ${esc(hit.resolution)} Å` : ''}<br><span style="color:var(--text-2)">${esc(hit.title || '')}</span></div><button class="btn btn-sm btn-primary" data-ligand-pdb-use="${i}">Load this PDB ID</button></div>`).join('');
    box.querySelector('[data-ligand-pdb-search]')?.addEventListener('click', () => search(statePdb.query));
    box.querySelectorAll('[data-ligand-pdb-use]').forEach(button => button.addEventListener('click', () => usePdb(Number(button.dataset.ligandPdbUse))));
  }

  function cancelSearch() {
    sequence++;
    if (activeSearch) activeSearch.abort();
    activeSearch = null;
    statePdb.loading = false;
    statePdb.error = 'Search cancelled.';
    render();
  }

  async function search(smiles) {
    const query = (smiles || '').trim(), seq = ++sequence;
    if (activeSearch) activeSearch.abort();
    activeSearch = null;
    statePdb.query = query;
    statePdb.hits = [];
    statePdb.error = '';
    if (!query) { statePdb.loading = false; render(); return; }
    statePdb.loading = true; render();
    const controller = new AbortController();
    activeSearch = controller;
    const timeout = setTimeout(() => controller.abort(), 24000);
    try {
      const response = await fetch(`/api/ligand-pdb/search?smiles=${encodeURIComponent(query)}&limit=6`, { signal: controller.signal, cache: 'no-store' });
      const data = await response.json().catch(() => ({}));
      if (!response.ok) throw new Error(data.detail || `Search failed (${response.status})`);
      if (seq !== sequence) return;
      statePdb.query = data.query_smiles || query;
      statePdb.hits = data.results || [];
    } catch (error) {
      if (seq !== sequence) return;
      statePdb.hits = [];
      statePdb.error = error.name === 'AbortError'
        ? 'RCSB search timed out. Press Find PDBs with this native ligand to retry.'
        : `RCSB search unavailable: ${error.message}`;
      statePdb.loading = false;
      render();
      return;
    } finally {
      clearTimeout(timeout);
      if (activeSearch === controller) activeSearch = null;
    }
    statePdb.loading = false;
    render();
  }

  async function usePdb(index) {
    const hit = statePdb.hits[index];
    if (!hit) return;
    try {
      setLoading(true, `Preparing PDB ${hit.pdb_id} and centring the docking box…`);
      const data = await apiFetch('/api/ligand-pdb/import', {
        method: 'POST',
        body: JSON.stringify({ pdb_id: hit.pdb_id, ccd_id: hit.ccd_id, smiles: statePdb.query })
      });
      const importedTarget = data.info || {};
      const matchedSmiles = statePdb.query;
      importedTarget.native_ligand = hit.ccd_id;
      importedTarget.native_ligand_removed = true;
      if (!Array.isArray(importedTarget.native_ligands) || !importedTarget.native_ligands.length) {
        // The ligand search already matched this CCD component by chemistry.
        // Carry that exact query into the imported target so native docking
        // does not depend on a second RCSB/CCD lookup succeeding.
        importedTarget.native_ligands = [{
          ccd_id: hit.ccd_id,
          name: hit.ligand_name || hit.ccd_id,
          smiles: matchedSmiles,
          atom_count: null,
          chain: '—',
          residue: '—',
          center: importedTarget.center,
          size: importedTarget.size,
        }];
      }
      state.targets[data.target_id] = importedTarget;
      state.pdbCache = state.pdbCache || {};
      delete state.pdbCache[data.target_id];
      state.selectedTarget = data.target_id;
      setSMILES(matchedSmiles);
      selectTarget(data.target_id);
      showToast(`Loaded PDB ${hit.pdb_id} with native ligand ${hit.ccd_id} removed.`, 'success', '🧬');
    } catch (error) {
      showToast(`Could not load PDB ${hit.pdb_id}: ${error.message}`, 'error', '❌');
    } finally {
      setLoading(false);
    }
  }

  function schedule(smiles) {
    clearTimeout(timer);
    const query = (smiles || '').trim();
    timer = setTimeout(() => {
      if (query === lastQuery) return;
      lastQuery = query;
      search(query);
    }, 650);
  }

  const originalSetSMILES = window.setSMILES;
  if (typeof originalSetSMILES === 'function') {
    window.setSMILES = function () {
      const result = originalSetSMILES.apply(this, arguments);
      schedule(typeof getSMILES === 'function' ? getSMILES() : arguments[0]);
      return result;
    };
  }
  document.addEventListener('input', event => {
    if (event.target && event.target.id === 'smiles-input') schedule(event.target.value);
  });
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', () => {
    render();
    const input = document.getElementById('smiles-input');
    if (input && input.value.trim()) schedule(input.value);
  });
  else {
    render();
    const input = document.getElementById('smiles-input');
    if (input && input.value.trim()) schedule(input.value);
  }
})();
// === END LIGAND_PDB_PICKER_V2 ===

