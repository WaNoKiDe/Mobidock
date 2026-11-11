// ══ CELL 46 · static/js/epoch.js ══
/* ─── epoch.js — epoch — competition epochs, WebSocket, credits, history, entropy engine ─── */


// ── Epoch System — WebSocket Client & UI Logic ───────────────────────────────

function initEpoch() {
  connectEpochWs();
  fetchEpochState();
  fetchPlayerCredits();
}

function connectEpochWs() {
  const protocol = location.protocol === 'https:' ? 'wss:' : 'ws:';
  const wsUrl = `${protocol}//${location.host}/ws/epoch`;
  const ws = new WebSocket(wsUrl);
  state.epochWs = ws;

  ws.onopen = () => {
    console.log('Epoch WS connected');
  };

  ws.onmessage = (evt) => {
    try {
      const msg = JSON.parse(evt.data);
      handleEpochMessage(msg);
    } catch (e) {
      console.warn('WS parse error:', e);
    }
  };

  ws.onclose = () => {
    console.warn('Epoch WS disconnected, reconnecting in 3s...');
    setTimeout(connectEpochWs, 3000);
  };

  ws.onerror = (err) => {
    console.error('Epoch WS error:', err);
  };

  // Ping keepalive
  setInterval(() => {
    if (ws.readyState === WebSocket.OPEN) {
      ws.send('ping');
    }
  }, 30000);
}

function handleEpochMessage(msg) {
  switch (msg.type) {
    case 'current_epoch':
    case 'new_epoch':
      state.epoch = msg.data?.epoch || msg.data;
      updateEpochUI();
      if (msg.type === 'new_epoch') {
        showToast(`New epoch! Target: ${state.epoch.protein_id} — Bounty: ${state.epoch.bounty_amount} CC`, 'info', '⚔️');
      }
      break;
    case 'epoch_tick':
      if (state.epoch && msg.data.epoch_id === state.epoch.epoch_id) {
        state.epoch.remaining_seconds = msg.data.remaining_seconds;
        state.epoch.solved_by = msg.data.solved_by;
        updateEpochTimer();
      }
      break;
    case 'epoch_solved':
      showEpochSolvedBanner(msg.data.solved_by);
      showToast(`Epoch ${msg.data.epoch_id} solved by ${msg.data.solved_by}!`, 'success', '🏆');
      break;
    case 'epoch_expired':
      showToast(`Epoch ${msg.data.epoch_id} expired.`, 'info', '⏰');
      break;
    case 'entropy_complete':
      handleEntropyComplete(msg.data);
      break;
    case 'pong':
      break;
    default:
      console.log('Unknown WS message:', msg);
  }
}

async function fetchEpochState() {
  try {
    const data = await apiFetch('/api/v1/epoch/target');
    if (data.active && data.epoch) {
      state.epoch = data.epoch;
      updateEpochUI();
    }
  } catch (e) {
    console.warn('Failed to fetch epoch:', e);
  }
}

function updateEpochUI() {
  const ep = state.epoch;
  if (!ep) return;

  document.getElementById('epoch-id').textContent = ep.epoch_id;
  document.getElementById('epoch-target').textContent = ep.protein_id;
  document.getElementById('epoch-ligand').textContent = `${ep.ligand_name} (${ep.ligand_smiles})`;
  document.getElementById('epoch-threshold').textContent = `${ep.difficulty_threshold_kcal_mol} kcal/mol`;
  document.getElementById('epoch-bounty').textContent = `${ep.bounty_amount} CC`;
  document.getElementById('epoch-status').textContent = ep.solved_by
    ? `✅ Solved by ${ep.solved_by}`
    : '⚡ Active — Solve Now!';

  updateEpochTimer();
  startEpochTimerInterval();
}

function updateEpochTimer() {
  const ep = state.epoch;
  if (!ep) return;

  const remaining = ep.remaining_seconds ?? 0;
  const total = ep.epoch_timer ?? 300;
  const fraction = Math.max(0, remaining / total);

  document.getElementById('timer-seconds').textContent = Math.ceil(remaining);

  // SVG circle animation (circumference = 2πr = 2π(54) ≈ 339.3)
  const circumference = 2 * Math.PI * 54;
  const fill = document.getElementById('timer-fill');
  if (fill) {
    fill.style.strokeDasharray = `${circumference}`;
    fill.style.strokeDashoffset = `${circumference * (1 - fraction)}`;
  }

  // Color shifts: green > yellow > red as time runs out
  const timerEl = document.getElementById('epoch-timer');
  if (fraction > 0.5) {
    timerEl.style.setProperty('--timer-color', '#00f5d4');
  } else if (fraction > 0.2) {
    timerEl.style.setProperty('--timer-color', '#ffd166');
  } else {
    timerEl.style.setProperty('--timer-color', '#ff6b6b');
  }
}

function startEpochTimerInterval() {
  if (state.epochTimerInterval) clearInterval(state.epochTimerInterval);
  state.epochTimerInterval = setInterval(() => {
    if (state.epoch && state.epoch.remaining_seconds > 0) {
      state.epoch.remaining_seconds = Math.max(0, state.epoch.remaining_seconds - 1);
      updateEpochTimer();
    }
  }, 1000);
}

function showEpochSolvedBanner(solver) {
  const banner = document.getElementById('epoch-solved-banner');
  document.getElementById('banner-solver').textContent = `by ${solver}`;
  banner.style.display = 'block';
  banner.classList.add('animate');
  setTimeout(() => {
    banner.style.display = 'none';
    banner.classList.remove('animate');
  }, 5000);
}

// ── Epoch Submission ──────────────────────────────────────────────────────

async function submitEpochSolution() {
  if (!state.epoch) {
    showToast('No active epoch', 'error', '⚠️');
    return;
  }

  const player = document.getElementById('player-input').value.trim() || 'Anonymous';
  const btn = document.getElementById('epoch-submit-btn');
  btn.disabled = true;
  btn.textContent = '⏳ Verifying...';

  try {
    const res = await apiFetch('/api/v1/epoch/submit', {
      method: 'POST',
      body: JSON.stringify({
        epoch_id: state.epoch.epoch_id,
        solver_id: player,
        client_timestamp: Date.now() / 1000,
        ligand_transform: {
          position: [0, 0, 0],
          rotation: [0, 0, 0, 1],
          torsion_angles: [],
        },
      }),
    });

    if (res.status === 'verified') {
      showToast(
        `✅ Verified! ${res.is_first_solve ? 'FIRST SOLVE! ' : ''}+${res.credits_awarded} CC | E=${res.verified_energy.toFixed(2)} kcal/mol`,
        'success', '🏆'
      );
      state.playerCredits = res.total_credits;
      state.lastPoseHash = res.pose_hash;
      state.lastEpochSolve = {
        epoch_id: state.epoch.epoch_id,
        target_id: state.epoch.protein_id,
        ligand_smiles: state.epoch.ligand_smiles,
        binding_energy: res.verified_energy,
        pose_hash: res.pose_hash,
      };
      updateCreditsDisplay();
      fetchPlayerHistory();
    } else {
      showToast(`❌ Rejected: ${res.reason}`, 'error', '⚠️');
    }
  } catch (err) {
    showToast(`Submission failed: ${err.message}`, 'error', '❌');
  } finally {
    btn.disabled = false;
    btn.textContent = '⚡ Submit Solution';
  }
}

// ── Credits Display ───────────────────────────────────────────────────────

async function fetchPlayerCredits() {
  const player = document.getElementById('player-input').value.trim() || 'Anonymous';
  try {
    const data = await apiFetch(`/api/v1/player/${encodeURIComponent(player)}/credits`);
    state.playerCredits = data.credits;
    updateCreditsDisplay();
  } catch (e) {
    console.warn('Failed to fetch credits:', e);
  }
}

function updateCreditsDisplay() {
  const headerCredits = document.getElementById('credits-amount');
  const dashCredits = document.getElementById('dashboard-credits');
  const amount = state.playerCredits;

  if (headerCredits) {
    const prev = parseInt(headerCredits.textContent) || 0;
    animateValue(headerCredits, prev, amount, 0, 600);
    // Pulse animation on credit gain
    if (amount > prev) {
      const badge = document.getElementById('credits-badge');
      badge.classList.add('pulse');
      setTimeout(() => badge.classList.remove('pulse'), 800);
    }
  }
  if (dashCredits) {
    const prev = parseInt(dashCredits.textContent) || 0;
    animateValue(dashCredits, prev, amount, 0, 800);
  }
}

// ── Player History ────────────────────────────────────────────────────────

async function fetchPlayerHistory() {
  const player = document.getElementById('player-input').value.trim() || 'Anonymous';
  try {
    const data = await apiFetch(`/api/v1/player/${encodeURIComponent(player)}/history`);
    renderResolutionLog(data.history || []);
    renderEntropySection();
  } catch (e) {
    console.warn('Failed to fetch history:', e);
  }
}

function renderResolutionLog(entries) {
  const container = document.getElementById('resolution-log');
  if (!entries.length) {
    container.innerHTML = '<div style="text-align:center; color:var(--text-3); font-size:12px; padding:20px">No epoch solves yet — compete to earn credits!</div>';
    return;
  }

  let html = '<div class="log-entries">';
  entries.forEach(e => {
    const date = new Date(e.timestamp * 1000).toLocaleString();
    const badge = e.is_first_solve ? '<span class="first-solve-badge">🥇 FIRST</span>' : '';
    const typeIcon = e.entry_type === 'entropy_bonus' ? '🔥' : e.entry_type === 'entropy_debit' ? '⚡' : '⚔️';
    const creditClass = e.credits_awarded >= 0 ? 'credit-positive' : 'credit-negative';
    html += `
      <div class="log-entry">
        <span class="log-icon">${typeIcon}</span>
        <div class="log-details">
          <div class="log-epoch">${e.epoch_id} — ${e.target_id} ${badge}</div>
          <div class="log-meta">${date} · E=${e.verified_binding_energy.toFixed(2)} kcal/mol</div>
        </div>
        <span class="log-credits ${creditClass}">${e.credits_awarded >= 0 ? '+' : ''}${e.credits_awarded} CC</span>
      </div>
    `;
  });
  html += '</div>';
  container.innerHTML = html;
}

// ── Entropy Engine UI ─────────────────────────────────────────────────────

function renderEntropySection() {
  const container = document.getElementById('entropy-jobs');
  if (!state.lastEpochSolve) {
    container.innerHTML = '<div style="text-align:center; color:var(--text-3); font-size:12px; padding:20px">Complete an epoch solve to unlock entropy testing.</div>';
    return;
  }

  let html = `
    <div class="entropy-submit-card">
      <div class="entropy-info">
        <strong>Last verified pose:</strong> ${state.lastEpochSolve.epoch_id} — ${state.lastEpochSolve.target_id}<br>
        <span class="font-mono" style="font-size:11px; color:var(--text-2)">E = ${state.lastEpochSolve.binding_energy.toFixed(2)} kcal/mol</span>
      </div>
      <button class="btn btn-primary" id="entropy-submit-btn" onclick="submitEntropyTest()" ${state.playerCredits < 100 ? 'disabled title="Need 100 CC"' : ''}>
        🔥 Run MD Test (100 CC)
      </button>
    </div>
  `;

  // Fetch and show existing jobs
  fetchEntropyJobs(html, container);
}

async function fetchEntropyJobs(prefixHtml, container) {
  const player = document.getElementById('player-input').value.trim() || 'Anonymous';
  try {
    const data = await apiFetch(`/api/v1/entropy/${encodeURIComponent(player)}/jobs`);
    let html = prefixHtml || '';
    if (data.jobs && data.jobs.length > 0) {
      data.jobs.forEach(job => {
        html += renderEntropyJobCard(job);
      });
    }
    container.innerHTML = html;
  } catch (e) {
    container.innerHTML = (prefixHtml || '') + '<div style="color:var(--text-3); font-size:12px; padding:8px">Could not load entropy jobs.</div>';
  }
}

function renderEntropyJobCard(job) {
  const statusIcons = {
    queued: '⏳', running: '⚙️', complete: '✅', failed: '❌'
  };
  const verdictClass = job.verdict === 'thermodynamically_verified' ? 'entropy-verified' :
                       job.verdict === 'false_equilibrium' ? 'entropy-rejected' : '';
  const verdictLabel = job.verdict === 'thermodynamically_verified' ? '✅ VERIFIED (3×)' :
                       job.verdict === 'false_equilibrium' ? '❌ FALSE EQUILIBRIUM' : '⏳ Pending';

  let rmsdBar = '';
  if (job.rmsd_trajectory && job.rmsd_trajectory.length > 0) {
    const maxRmsd = Math.max(...job.rmsd_trajectory, 3);
    const barPoints = job.rmsd_trajectory.map((r, i) => {
      const x = (i / (job.rmsd_trajectory.length - 1)) * 100;
      const y = 100 - (r / maxRmsd) * 100;
      return `${x},${y}`;
    }).join(' ');
    rmsdBar = `
      <div class="rmsd-chart">
        <svg viewBox="0 0 100 100" preserveAspectRatio="none">
          <line x1="0" y1="${100 - (2.0/maxRmsd)*100}" x2="100" y2="${100 - (2.0/maxRmsd)*100}" stroke="#ff6b6b" stroke-width="0.5" stroke-dasharray="2,2"/>
          <polyline points="${barPoints}" fill="none" stroke="${job.verdict === 'thermodynamically_verified' ? '#00f5d4' : '#ff6b6b'}" stroke-width="1.5"/>
        </svg>
        <div class="rmsd-labels">
          <span>0 ns</span>
          <span style="color:var(--text-3)">RMSD threshold: 2.0 Å</span>
          <span>${job.rmsd_trajectory.length} ns</span>
        </div>
      </div>
    `;
  }

  return `
    <div class="entropy-card ${verdictClass}">
      <div class="entropy-card-header">
        <span class="entropy-job-id">${statusIcons[job.status] || '⚙️'} ${job.job_id}</span>
        <span class="entropy-verdict-badge ${verdictClass}">${verdictLabel}</span>
      </div>
      <div class="entropy-card-body">
        <div class="entropy-stat">Mean RMSD: <strong>${job.mean_rmsd.toFixed(2)} Å</strong></div>
        <div class="entropy-stat">Max RMSD: <strong>${job.max_rmsd.toFixed(2)} Å</strong></div>
        ${job.credits_bonus > 0 ? `<div class="entropy-stat credit-positive">Bonus: +${job.credits_bonus} CC</div>` : ''}
      </div>
      ${rmsdBar}
    </div>
  `;
}

async function submitEntropyTest() {
  if (!state.lastEpochSolve) {
    showToast('Complete an epoch solve first', 'error', '⚠️');
    return;
  }
  if (state.playerCredits < 100) {
    showToast('Insufficient credits (need 100 CC)', 'error', '⚠️');
    return;
  }

  const player = document.getElementById('player-input').value.trim() || 'Anonymous';
  const btn = document.getElementById('entropy-submit-btn');
  if (btn) { btn.disabled = true; btn.textContent = '⏳ Submitting...'; }

  try {
    const res = await apiFetch('/api/v1/entropy/submit', {
      method: 'POST',
      body: JSON.stringify({
        user_id: player,
        epoch_id: state.lastEpochSolve.epoch_id,
        target_id: state.lastEpochSolve.target_id,
        pose_hash: state.lastEpochSolve.pose_hash,
        ligand_smiles: state.lastEpochSolve.ligand_smiles,
        binding_energy: state.lastEpochSolve.binding_energy,
      }),
    });

    showToast(`MD job ${res.job_id} queued! -${res.credits_deducted} CC`, 'info', '🔥');
    state.playerCredits = res.remaining_credits;
    updateCreditsDisplay();

    // Poll for completion
    pollEntropyJob(res.job_id);
  } catch (err) {
    showToast(`Entropy submit failed: ${err.message}`, 'error', '❌');
  } finally {
    if (btn) { btn.disabled = false; btn.textContent = '🔥 Run MD Test (100 CC)'; }
  }
}

async function pollEntropyJob(jobId) {
  let attempts = 0;
  const poll = setInterval(async () => {
    attempts++;
    try {
      const job = await apiFetch(`/api/v1/entropy/job/${jobId}`);
      if (job.status === 'complete' || job.status === 'failed') {
        clearInterval(poll);
        handleEntropyComplete(job);
      }
    } catch (e) {
      if (attempts > 30) clearInterval(poll);
    }
  }, 2000);
}

function handleEntropyComplete(job) {
  if (job.verdict === 'thermodynamically_verified') {
    showToast(`✅ THERMODYNAMICALLY VERIFIED! +${job.credits_bonus} CC bonus (3× multiplier)`, 'success', '🔥');
  } else if (job.verdict === 'false_equilibrium') {
    showToast(`❌ False Equilibrium — pose unstable (RMSD ${job.mean_rmsd.toFixed(2)}Å)`, 'error', '🔥');
  }
  fetchPlayerCredits();
  renderEntropySection();
}
