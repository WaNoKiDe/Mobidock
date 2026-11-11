// ══ CELL 45 · static/js/system.js ══
/* ─── system.js — system — JSME molecule editor setup and server health check ─── */


// ── JSME setup ────────────────────────────────────────────────────────────────
function initJSME() {
  let built = false, tries = 0;
  function build() {
    if (built) return true;
    if (!(window.JSApplet && window.JSApplet.JSME)) return false;
    try {
      window.jsmeApp = new JSApplet.JSME('jsme-container', '100%', '320px', { options: 'query,hydrogens' });
      window.jsmeApp.setCallBack('AfterStructureModified', function () {
        if (state.jsmeSyncLock) return;
        try {
          const s = window.jsmeApp.smiles();
          if (s && s.trim() && s !== '.') document.getElementById('smiles-input').value = s.trim();
        } catch (e) {}
      });
      built = true;
      setTimeout(() => setSMILES(getSMILES() || 'CC(=O)Oc1ccccc1C(=O)O'), 500);
    } catch (e) { console.warn('JSME init failed:', e); }
    return built;
  }
  window.jsmeOnLoad = build;                 // JSME calls this when it finishes loading
  const timer = setInterval(() => {          // ...and we also poll, in case it loaded first
    tries++;
    if (build()) { clearInterval(timer); return; }
    if (tries === 30) showToast('2D editor still loading (needs access to jsme-editor.github.io). You can type a SMILES meanwhile.', 'warning', '⚠️');
    if (tries > 240) clearInterval(timer);
  }, 500);
}

// ── Health check ──────────────────────────────────────────────────────────────
async function checkHealth() {
  try {
    const h = await apiFetch('/api/health');
    const dot = document.getElementById('status-dot');
    const txt = document.getElementById('status-text');
    if (h.status === 'ok') {
      dot.style.background = 'var(--green)';
      const ready = Object.values(h.engines || {}).filter(v => v === 'available').length;
      const total = Object.keys(h.engines || {}).length;
      const modeLabel = h.docking_mode === 'none' ? '⚠ Vina missing' : '⚗️ Vina';
      txt.textContent = `${modeLabel} | ${ready}/${total} engines ready${h.gpu ? ' | GPU' : ''} | RDKit ${h.rdkit}`;
    }
  } catch(e) {
    document.getElementById('status-dot').style.background = 'var(--red)';
    document.getElementById('status-text').textContent = 'Server offline';
  }
}
