// ══ CELL 36 · static/js/engines.js ══
/* ─── engines.js — engine picker + setup panel ────────────────────────────────
 * Every engine is real. If one is not ready, clicking it opens a panel that
 * lists exactly what is missing (and, for DiffDock-L, lets you paste your free
 * NVIDIA API key right here).
 */

async function loadEngines() {
  try {
    state.engines = await apiFetch('/api/engines');
  } catch (e) {
    console.warn('Failed to load engines:', e);
    state.engines = [];
  }
  renderEngineSelector();
}

function getDiffDockBrowserKey() {
  try {
    const key = window.localStorage.getItem('mobidock.diffdock.nim-key') || '';
    return key.startsWith('nvapi-') && key.length >= 20 ? key : '';
  } catch (_) { return ''; }
}

function engineRequiresNvidiaGpu(eng) {
  return eng?.name === 'diffdock' && (eng.problems || []).some(problem => problem.includes('only on NVIDIA GPU laptops'));
}

function renderEngineSelector() {
  const wrap = document.getElementById('engine-options');
  wrap.innerHTML = '';
  state.engineChoice = state.engineChoice || state.selectedEngine || 'vina';
  const all = document.createElement('div');
  all.className = `engine-option${state.engineChoice === 'all' ? ' active' : ''}`;
  all.dataset.engine = 'all';
  const cpuOnlyLaptop = (state.engines || []).some(engineRequiresNvidiaGpu);
  all.innerHTML = `<span class="engine-icon">⚗️</span><span class="engine-name">${cpuOnlyLaptop ? 'All three CPU engines' : 'All four engines'}</span><span class="engine-status available">RUN ALL</span>`;
  all.title = (state.engines || []).map(engine => `${engine.display_name}: ${engine.status === 'available' ? 'ready' : (engine.problems || []).join('; ')}`).join('\n');
  all.addEventListener('click', () => {
    state.engineChoice = 'all';
    document.querySelectorAll('.engine-option').forEach(option => option.classList.toggle('active', option.dataset.engine === 'all'));
    const setup = document.getElementById('engine-setup');
    if (setup) setup.style.display = 'none';
  });
  wrap.appendChild(all);
  state.engines.forEach(eng => {
    const blocked = engineRequiresNvidiaGpu(eng);
    const ready = !blocked && (eng.status === 'available' || (eng.name === 'diffdock' && !!getDiffDockBrowserKey()));
    const opt = document.createElement('div');
    opt.className = `engine-option${state.engineChoice === eng.name ? ' active' : ''}`;
    opt.dataset.engine = eng.name;
    opt.style.setProperty('--eng-color', eng.color);
    const statusLabel = blocked ? 'GPU REQUIRED' : (eng.status === 'available' ? '✓ READY' : (ready ? '✓ MY KEY' : '⚙ SETUP'));
    opt.innerHTML = `
      <span class="engine-icon">${eng.icon}</span>
      <span class="engine-name">${eng.display_name}</span>
      <span class="engine-status ${ready ? 'available' : 'needs_config'}">${statusLabel}</span>`;
    opt.title = eng.description + (ready ? '' : '\n\nMissing:\n- ' + eng.problems.join('\n- '));
    opt.addEventListener('click', () => {
      state.selectedEngine = eng.name;
      state.engineChoice = eng.name;
      document.querySelectorAll('.engine-option').forEach(o => o.classList.toggle('active', o === opt));
      renderEngineSetupPanel(eng);
    });
    wrap.appendChild(opt);
  });
  const current = state.engineChoice === 'all' ? null : state.engines.find(e => e.name === state.engineChoice);
  if (current) renderEngineSetupPanel(current);
  else {
    const setup = document.getElementById('engine-setup');
    if (setup) setup.style.display = 'none';
  }
}

function renderEngineSetupPanel(eng) {
  let panel = document.getElementById('engine-setup');
  if (!panel) {
    panel = document.createElement('div');
    panel.id = 'engine-setup';
    document.getElementById('engine-options').insertAdjacentElement('afterend', panel);
  }
  if (eng.status === 'available') { panel.style.display = 'none'; panel.innerHTML = ''; return; }
  const userKeySaved = eng.name === 'diffdock' && !!getDiffDockBrowserKey();
  const gpuBlocked = engineRequiresNvidiaGpu(eng);

  const keyBox = eng.name === 'diffdock' && !gpuBlocked ? `
    <div style="display:flex;gap:6px;margin-top:8px">
      <input id="nim-key-input" type="password" class="smiles-input" placeholder="nvapi-…  (remembered in this browser)"
             style="flex:1;padding:6px;font-size:12px" autocomplete="off">
      <button class="btn btn-sm btn-primary" onclick="saveDiffDockKey()">Save key</button>
    </div>` : '';

  panel.style.cssText = 'display:block;margin-top:10px;padding:10px 12px;border-radius:8px;font-size:12px;' +
    'background:var(--bg-1);border:1px solid var(--border);color:var(--text-2)';
  panel.innerHTML = `
    <div style="color:var(--accent-4);font-weight:600;margin-bottom:4px">${eng.icon} ${eng.display_name} ${gpuBlocked ? 'requires an NVIDIA GPU' : (userKeySaved ? 'key is saved in this browser' : 'needs setup')}</div>
    ${gpuBlocked ? `<ul style="margin:0 0 4px 16px;padding:0">${eng.problems.map(p => `<li>${escapeHtml(p)}</li>`).join('')}</ul>` : (userKeySaved ? '<div>Your key stays on this device and is sent only with your DiffDock requests. Other visitors cannot use it.</div><button class="btn btn-sm" style="margin-top:6px" onclick="clearDiffDockKey()">Clear key</button>' : `<ul style="margin:0 0 4px 16px;padding:0">${eng.problems.map(p => `<li>${escapeHtml(p)}</li>`).join('')}</ul>`)}
    <div style="color:var(--text-3)">${escapeHtml(eng.setup_hint || '')}</div>
    ${keyBox}
    <button class="btn btn-sm" style="margin-top:8px" onclick="loadEngines()">↻ Re-check</button>`;
}

async function saveDiffDockKey() {
  const key = document.getElementById('nim-key-input').value.trim();
  if (!key) { showToast('Paste your NVIDIA API key first', 'error', '⚠️'); return; }
  if (!key.startsWith('nvapi-') || key.length < 20) { showToast('That does not look like an NVIDIA API key', 'error', '⚠️'); return; }
  try { window.localStorage.setItem('mobidock.diffdock.nim-key', key); }
  catch (_) { showToast('This browser could not save the key on this device', 'error', '❌'); return; }
  showToast('DiffDock key saved in this browser', 'success', '🧬');
  await loadEngines();
}

async function clearDiffDockKey() {
  try { window.localStorage.removeItem('mobidock.diffdock.nim-key'); } catch (_) {}
  await loadEngines();
}
