// ══ CELL 37 · static/js/gridbox.js ══
/* ─── gridbox.js — grid box inputs (centre X/Y/Z + size X/Y/Z) ────────────────
 * The box is a cuboid: three independent edge lengths in Å. Changes are sent to
 * the server (so every engine docks inside it) and redrawn in the 3D viewer.
 */

function readBoxInputs() {
  const num = (id, d) => { const v = parseFloat(document.getElementById(id).value); return isFinite(v) ? v : d; };
  return {
    center: [num('box-cx', 0), num('box-cy', 0), num('box-cz', 0)],
    size:   [Math.max(4, num('box-sx', 20)), Math.max(4, num('box-sy', 20)), Math.max(4, num('box-sz', 20))],
  };
}

function fillBoxInputs(center, size) {
  const set = (id, v, d) => { document.getElementById(id).value = parseFloat(v).toFixed(d); };
  set('box-cx', center[0], 2); set('box-cy', center[1], 2); set('box-cz', center[2], 2);
  set('box-sx', size[0], 1);   set('box-sy', size[1] ?? size[0], 1); set('box-sz', size[2] ?? size[0], 1);
}

function updateBoxDebounced() {
  clearTimeout(state.boxTimer);
  state.boxTimer = setTimeout(updateBox, 350);
}

function setCubicBox() {
  const x = document.getElementById('box-sx').value;
  document.getElementById('box-sy').value = x;
  document.getElementById('box-sz').value = x;
  updateBox();
}

async function updateBox() {
  if (!state.selectedTarget) return;
  const { center, size } = readBoxInputs();
  try {
    await apiFetch(`/api/target/${state.selectedTarget}/box`, {
      method: 'POST',
      body: JSON.stringify({ center, size }),
    });
    const t = state.targets[state.selectedTarget];
    if (t) { t.center = center; t.size = size; }
    if (state.viewer) {
      drawGridBox(state.viewer, center, size, !!state.ligModel);
      state.viewer.render();
    }
  } catch (e) {
    showToast(`Failed to update grid box: ${e.message}`, 'error', '❌');
  }
}
