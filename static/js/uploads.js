// ══ CELL 38 · static/js/uploads.js ══
/* ─── uploads.js — uploads — ligand/receptor file upload, fetching the target list ─── */


// ── Upload Handlers ───────────────────────────────────────────────────────────
function handleLigandUpload(event) {
  const file = event.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = function(e) {
    const content = e.target.result;
    if (window.jsmeApp && typeof window.jsmeApp.readMolFile === 'function') {
      state.jsmeSyncLock = true;
      window.jsmeApp.readMolFile(content);
      setTimeout(() => { state.jsmeSyncLock = false; }, 300);
      const s = window.jsmeApp.smiles();
      if (s) {
        document.getElementById('smiles-input').value = s;
        showToast('Ligand loaded successfully', 'success', '🧪');
      } else {
        showToast('Could not parse ligand', 'error', '⚠️');
      }
    } else {
      showToast('JSME editor not ready', 'error', '⚠️');
    }
  };
  reader.readAsText(file);
}

async function handleReceptorUpload(event) {
  const file = event.target.files[0];
  if (!file) return;

  const formData = new FormData();
  formData.append("file", file);

  setLoading(true, 'Uploading and processing custom PDB...');
  try {
    const res = await fetch('/api/upload-receptor', { method: 'POST', body: formData });
    const data = await res.json();
    if (!res.ok) throw new Error(data.detail || 'Upload failed');

    showToast('Receptor processed successfully!', 'success', '✅');
    await fetchTargets();
    selectTarget(data.target_id);
  } catch (err) {
    showToast(`Upload failed: ${err.message}`, 'error', '❌');
  } finally {
    setLoading(false);
  }
}

async function fetchTargets() {
  try {
    state.targets = await apiFetch('/api/targets');
    renderTargets();
    renderLeaderboardTabs();
  } catch(e) {
    console.warn('Failed to fetch targets:', e);
  }
}
