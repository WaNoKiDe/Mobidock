// ══ CELL 34 · static/js/core.js ══
/* ─── MedChem Game v3 — Frontend Application Logic ───────────────────────────────── */

// ── State ─────────────────────────────────────────────────────────────────────
const state = {
  selectedTarget: null,  // selected from an RCSB ligand-PDB result
  targets: {},
  lastResult: null,
  viewer: null,
  proteinModel: null,
  ligModel: null,
  boxModel: null,
  lbCurrentTab: null,
  selectedEngine: 'vina',
  engineChoice: 'vina',
  engines: [],
  dockingRuns: [],
  dockingRunSequence: 0,
  dockingBatchSequence: 0,
  dockingProgress: '',
  dockingActive: false,
  dockingPauseRequested: false,
  dockingStopRequested: false,
  dockingPaused: false,
  dockingStopped: false,
  dockingResumeBatch: null,
  dockingLoaderHidden: false,
  dockingCurrentEngine: '',
  dockingActiveContext: null,
  loadingKind: '',
  dockingError: '',
  dockingRmsdBest: {},
  dockingRmsdEpoch: {},
  pdbSearchPage: 1,
  pdbSearchQuery: '',
  pdbSearchTotal: 0,
  hoverSphere: null,
  boxShapes: [],
  interactionShapes: [],
  poseInteractionData: null,
  poseInteractionError: '',
  poseInteractionsLoading: false,
  interactionsEnabled: true,
  spinning: false,
  boxTimer: null,
  // v3 epoch state
  epoch: null,
  epochWs: null,
  epochTimerInterval: null,
  playerCredits: 0,
  lastPoseHash: null,
  lastEpochSolve: null,
};

// ── Quick-fill molecule library ───────────────────────────────────────────────
const QUICK_MOLS = [
  { name: 'Aspirin',     smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'Ibuprofen',   smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
  { name: 'Caffeine',    smiles: 'Cn1c(=O)c2c(ncn2C)n(c1=O)C' },
  { name: 'Benzene',     smiles: 'c1ccccc1' },
  { name: 'Dopamine',    smiles: 'NCCc1ccc(O)c(O)c1' },
  { name: 'Serotonin',   smiles: 'NCCc1c[nH]c2ccc(O)cc12' },
  { name: 'Melatonin',   smiles: 'COc1ccc2[nH]cc(CCNC(C)=O)c2c1' },
  { name: 'Resveratrol', smiles: 'Oc1ccc(cc1)/C=C/c1cc(O)cc(O)c1' },
];

// ── API helpers ───────────────────────────────────────────────────────────────
async function apiFetch(path, opts = {}) {
  const res = await fetch(path, { headers: { 'Content-Type': 'application/json' }, ...opts });
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: res.statusText }));
    throw new Error(err.detail || 'API error');
  }
  return res.json();
}

// ── Toast notifications & Notification Center ──────────────────────────────────
let notificationHistory = [];

function updateNotificationDropdown() {
  const badge = document.getElementById('notif-badge');
  const list = document.getElementById('notif-list');
  if (!badge || !list) return;

  if (notificationHistory.length > 0) {
    // Only show badge if the dropdown is NOT currently open
    const dropdown = document.getElementById('notif-dropdown');
    if (dropdown && dropdown.style.display !== 'block') {
      badge.style.display = 'block';
      badge.innerText = notificationHistory.length > 99 ? '99+' : notificationHistory.length;
    }

    list.innerHTML = notificationHistory.map(n => `
      <div style="background:var(--surface-1); padding:10px; border-radius:6px; font-size:0.85rem; border-left:3px solid ${n.type==='error'?'var(--red)':n.type==='success'?'var(--green)':'var(--accent-1)'};">
        <div style="display:flex; justify-content:space-between; margin-bottom:4px;">
          <span style="font-size:1.1rem;">${n.icon}</span>
          <span style="color:var(--text-3); font-size:0.75rem;">${n.time.toLocaleTimeString()}</span>
        </div>
        <div style="color:var(--text-1); line-height:1.4;">${n.msg}</div>
      </div>
    `).join('');
  } else {
    list.innerHTML = '<div style="text-align:center; color:var(--text-3); font-size:0.9rem; padding:10px;">No notifications yet.</div>';
  }
}

function showToast(msg, type = 'info', icon = '🔬') {
  // Push to history array (newest first)
  notificationHistory.unshift({msg, type, icon, time: new Date()});
  updateNotificationDropdown();

  const container = document.getElementById('toast-container');
  // Clear any existing toasts to ensure ONLY 1 is shown at a time
  container.innerHTML = '';

  const toast = document.createElement('div');
  toast.className = `toast ${type}`;
  toast.innerHTML = `<span style="font-size:18px">${icon}</span> <span>${msg}</span>`;
  container.appendChild(toast);

  setTimeout(() => {
    toast.style.animation = 'slide-in .3s ease reverse';
    setTimeout(() => { if (toast.parentNode) toast.remove(); }, 300);
  }, 3500);
}

// ── Loading overlay ───────────────────────────────────────────────────────────
function setLoading(show, step = '', kind = 'general') {
  const overlay = document.getElementById('loading-overlay');
  if (!overlay) return;
  state.loadingKind = show ? kind : '';
  document.body.classList.toggle('non-docking-loading', Boolean(show && kind !== 'docking'));
  const isDockingLoad = show && kind === 'docking' && state.dockingActive;
  overlay.classList.toggle('show', Boolean(show && !(isDockingLoad && state.dockingLoaderHidden)));
  const dockingControls = document.getElementById('docking-loading-controls');
  if (dockingControls) dockingControls.hidden = !isDockingLoad;
  if (step) {
    const label = document.getElementById('loading-step');
    if (label) label.textContent = step;
  }
}

// ── Animated counter ──────────────────────────────────────────────────────────
function animateValue(el, from, to, decimals = 2, duration = 800) {
  const start = performance.now();
  const update = (now) => {
    const t = Math.min((now - start) / duration, 1);
    const ease = t < 0.5 ? 2*t*t : -1+(4-2*t)*t;
    const val = from + (to - from) * ease;
    el.textContent = val.toFixed(decimals);
    if (t < 1) requestAnimationFrame(update);
  };
  requestAnimationFrame(update);
}

// ── Tab Navigation ────────────────────────────────────────────────────────────
function initTabs() {
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const tabId = btn.dataset.tab;
      // Update buttons
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      // Update panels
      document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
      const panel = document.getElementById(`panel-${tabId}`);
      if (panel) panel.classList.add('active');
    });
  });
}

// ── Safe HTML text ────────────────────────────────────────────────────────────
function escapeHtml(str) {
  return String(str).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;').replace(/'/g,'&#39;');
}
