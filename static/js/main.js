// ══ CELL 47 · static/js/main.js ══
/* ─── main.js — main — init() and page start-up (loaded LAST) ─── */


// ── Main init ─────────────────────────────────────────────────────────────────
async function init() {
  // Keep the Dock control responsive even when the server check fails early.
  // The click then reports the connection problem instead of remaining inert.
  const dockButton = document.getElementById('dock-btn');
  if (dockButton && !dockButton.dataset.dockHandlerBound) {
    dockButton.dataset.dockHandlerBound = '1';
    dockButton.addEventListener('click', async () => {
      if (!state.targets || !state.targets[state.selectedTarget]) {
        showToast('MOBIDock is not connected to a loaded target. Restart the updated notebook server and open its link.', 'error', '🔴');
        return;
      }
      try { await runDocking(); }
      catch (error) {
        console.error('Docking could not start:', error);
        showToast(`Docking could not start: ${error.message || error}`, 'error', '❌');
      }
    });
  }
  const smilesInput = document.getElementById('smiles-input');
  if (smilesInput && !smilesInput.dataset.dockHandlerBound) {
    smilesInput.dataset.dockHandlerBound = '1';
    smilesInput.addEventListener('keydown', e => {
      if (e.key === 'Enter') runDocking();
    });
  }

  try { initBackground(); } catch (e) { console.warn('background failed', e); }
  initJSME();
  initTabs();
  renderQuickFills();

  try {
    const schema = await apiFetch('/openapi.json');
    const paths = schema.paths || {};
    const required = ['/api/dock', '/api/dock/rmsd', '/api/targets/{target_id}/native-ligands', '/api/ligand-pdb/import', '/api/engines', '/api/synth/explain3d', '/api/docking-results/export.xlsx'];
    const missing = required.filter(path => !paths[path]);
    if (missing.length) {
      state.serverCompatibilityWarning = `This tab is connected to an older MobiDock server (missing ${missing.join(', ')}). Restart from the updated MobiDock_Laptop notebook and open the link it displays.`;
      renderDockingRuns();
      showToast('This MobiDock tab is using an older server. Open the link from the updated notebook.', 'error', '⚠️');
    } else state.serverCompatibilityWarning = '';
  } catch (error) {
    state.serverCompatibilityWarning = `Could not verify the server API: ${error.message}`;
    showToast(`Could not verify the MobiDock server: ${error.message}`, 'error', '🔴');
  }

  // Load targets
  try {
    state.targets = await apiFetch('/api/targets');
    const targetIds = Object.keys(state.targets);
    renderTargets();

    // Never open a protein implicitly. The user selects one from the loaded-PDB list.
    if (state.selectedTarget && state.targets[state.selectedTarget]) {
      state.lbCurrentTab = state.selectedTarget;
      renderLeaderboardTabs();
      refreshLeaderboard();
      selectTarget(state.selectedTarget);
    } else {
      state.selectedTarget = null;
      state.lbCurrentTab = null;
      clearProteinViewer();
      renderTargetStructureDetails(null);
    }
  } catch(e) {
    showToast('Could not reach backend. Is the server running?', 'error', '🔴');
  }

  // Load engines
  await loadEngines();

  // Wire events
  const notifCenter = document.getElementById('notif-center');
  const notifDropdown = document.getElementById('notif-dropdown');
  const notifBadge = document.getElementById('notif-badge');
  const notifClearBtn = document.getElementById('notif-clear-btn');

  if (notifCenter && notifDropdown) {
    notifCenter.addEventListener('click', (e) => {
      e.stopPropagation();
      const isVisible = notifDropdown.style.display === 'block';
      notifDropdown.style.display = isVisible ? 'none' : 'block';
      if (!isVisible && notifBadge) {
        notifBadge.style.display = 'none'; // clear badge when opened
      }
    });
    document.addEventListener('click', (e) => {
      if (!notifCenter.contains(e.target)) notifDropdown.style.display = 'none';
    });
    if (notifClearBtn) {
      notifClearBtn.addEventListener('click', (e) => {
        e.stopPropagation();
        notificationHistory = [];
        updateNotificationDropdown();
        if (notifBadge) notifBadge.style.display = 'none';
      });
    }
  }

  // PDB search enter key
  document.getElementById('pdb-search-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') searchPDB();
  });
  // Ligand search enter key
  document.getElementById('ligand-search-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') searchLigands();
  });

  checkHealth();

  // v3: Initialize epoch system
  initEpoch();

  // Poll leaderboard every 15s
  setInterval(refreshLeaderboard, 15000);
}

document.addEventListener('DOMContentLoaded', init);

