// ══ CELL 40 · static/js/browsers.js ══
/* ─── browsers.js — browsers — PubChem search, quick-fills, protein browser (RCSB), ligand browser ─── */


// ── PubChem Quick Search (legacy inline) ──────────────────────────────────────
async function searchPubChem() {
  const input = document.getElementById('pubchem-search-input');
  const name = input.value.trim();
  if (!name) { showToast('Enter a chemical name to search', 'warning', '⚠️'); return; }

  const btn = document.getElementById('pubchem-search-btn');
  const originalText = btn.innerHTML;
  btn.innerHTML = '⏳...';
  btn.disabled = true;

  try {
    const res = await fetch(`/api/pubchem/search?name=${encodeURIComponent(name)}`);
    const data = await res.json();
    if (!res.ok) throw new Error(data.detail || 'PubChem search failed');

    setSMILES(data.smiles);
    showToast(`Found ${data.name}! SMILES loaded.`, 'success', '✅');
  } catch (err) {
    showToast(err.message, 'error', '❌');
  } finally {
    btn.innerHTML = originalText;
    btn.disabled = false;
  }
}

// ── PDB Import ────────────────────────────────────────────────────────────────
async function importPdb(pdbId) {
  if (!pdbId) {
    const input = document.getElementById('pdb-import-input');
    pdbId = input.value.trim();
  }
  if (!pdbId || pdbId.length !== 4) {
    showToast('Enter a valid 4-character PDB ID', 'warning', '⚠️');
    return;
  }

  setLoading(true, `Importing ${pdbId.toUpperCase()} from RCSB PDB...`);
  try {
    const res = await fetch('/api/pdb/import', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ pdb_id: pdbId })
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.detail || 'Import failed');

    await fetchTargets();
    selectTarget(data.target_id);
    const src = data.info && data.info.box_source ? ` Grid box: ${data.info.box_source}.` : '';
    showToast(`Target ${pdbId.toUpperCase()} imported and selected.${src}`, 'success', '✅');

    // Switch to dock tab
    document.querySelector('[data-tab="dock"]').click();
    return data;
  } catch (err) {
    showToast(`Import failed: ${err.message}`, 'error', '❌');
  } finally {
    setLoading(false);
  }
}


// ── Render quick-fill chips ───────────────────────────────────────────────────
function renderQuickFills() {
  const wrap = document.getElementById('quick-fills');
  wrap.innerHTML = '';
  QUICK_MOLS.forEach(m => {
    const chip = document.createElement('span');
    chip.className = 'qf-chip';
    chip.textContent = m.name;
    chip.setAttribute('title', m.smiles);
    chip.addEventListener('click', () => {
      setSMILES(m.smiles);
      showToast(`Loaded ${m.name}`, 'info', '🧪');
    });
    wrap.appendChild(chip);
  });
}

// ── Protein Browser (RCSB PDB) ────────────────────────────────────────────────
async function searchPDB() {
  const input = document.getElementById('pdb-search-input');
  const query = input.value.trim();
  if (!query) { showToast('Enter a search term', 'warning', '⚠️'); return; }

  state.pdbSearchQuery = query;
  state.pdbSearchPage = 1;
  await performPDBSearch();
}

async function performPDBSearch() {
  const results = document.getElementById('pdb-results');
  const pagination = document.getElementById('pdb-pagination');

  results.innerHTML = '<div class="browser-placeholder"><div class="spinner" style="width:40px;height:40px"></div></div>';

  try {
    const data = await apiFetch(`/api/pdb/search?q=${encodeURIComponent(state.pdbSearchQuery)}&page=${state.pdbSearchPage}&page_size=8`);

    state.pdbSearchTotal = data.total || 0;
    results.innerHTML = '';

    if (!data.results || data.results.length === 0) {
      results.innerHTML = '<div class="browser-placeholder"><div style="color:var(--text-3)">No results found. Try a different search term.</div></div>';
      pagination.style.display = 'none';
      return;
    }

    data.results.forEach(r => {
      const card = document.createElement('div');
      card.className = 'browser-result-card';
      card.innerHTML = `
        <div class="result-icon">🧬</div>
        <div class="result-info">
          <div class="result-title">${escapeHtml(r.title || 'Unknown Structure')}</div>
          <div class="result-subtitle">${escapeHtml(r.organism || 'Unknown organism')}</div>
          <div class="result-meta">
            <span class="result-tag">${r.pdb_id}</span>
            ${r.resolution ? `<span class="result-tag">${r.resolution}Å</span>` : ''}
            ${r.method ? `<span class="result-tag">${r.method}</span>` : ''}
          </div>
        </div>
        <button class="result-action">⬇ Import</button>
      `;
      card.querySelector('.result-action').addEventListener('click', (ev) => { ev.stopPropagation(); importPdb(r.pdb_id); });
      card.addEventListener('click', () => importPdb(r.pdb_id));
      results.appendChild(card);
    });

    // Pagination
    const totalPages = Math.ceil(state.pdbSearchTotal / 8);
    if (totalPages > 1) {
      pagination.style.display = 'flex';
      document.getElementById('pdb-page-info').textContent = `Page ${state.pdbSearchPage} of ${totalPages} (${state.pdbSearchTotal} results)`;
      document.getElementById('pdb-prev-btn').disabled = state.pdbSearchPage <= 1;
      document.getElementById('pdb-next-btn').disabled = state.pdbSearchPage >= totalPages;
    } else {
      pagination.style.display = 'none';
    }

  } catch (err) {
    results.innerHTML = `<div class="browser-placeholder"><div style="color:var(--red)">Search failed: ${err.message}</div></div>`;
  }
}

function pdbChangePage(delta) {
  state.pdbSearchPage += delta;
  performPDBSearch();
}

// ── Ligand Browser (PubChem) ──────────────────────────────────────────────────
async function searchLigands() {
  const input = document.getElementById('ligand-search-input');
  const query = input.value.trim();
  if (!query) { showToast('Enter a compound name', 'warning', '⚠️'); return; }

  const results = document.getElementById('ligand-results');
  results.innerHTML = '<div class="browser-placeholder"><div class="spinner" style="width:40px;height:40px"></div></div>';

  try {
    const data = await apiFetch(`/api/ligand/search?q=${encodeURIComponent(query)}`);

    results.innerHTML = '';

    if (!data.results || data.results.length === 0) {
      results.innerHTML = '<div class="browser-placeholder"><div style="color:var(--text-3)">No compounds found.</div></div>';
      return;
    }

    data.results.forEach(r => {
      const card = document.createElement('div');
      card.className = 'browser-result-card';
      card.innerHTML = `
        ${r.image_url ? `<img class="result-thumb" src="${r.image_url}" alt="${escapeHtml(r.name)}" onerror="this.style.display='none'">` : '<div class="result-icon">💊</div>'}
        <div class="result-info">
          <div class="result-title">${escapeHtml(r.name)}</div>
          <div class="result-subtitle font-mono" style="font-size:10px;max-width:300px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap">${escapeHtml(r.smiles || '')}</div>
          <div class="result-meta">
            ${r.cid ? `<span class="result-tag">CID: ${r.cid}</span>` : ''}
            ${r.formula ? `<span class="result-tag">${r.formula}</span>` : ''}
            ${r.mw ? `<span class="result-tag">${r.mw} Da</span>` : ''}
          </div>
        </div>
        <button class="result-action">⚡ Load</button>
      `;
      card.querySelector('.result-action').addEventListener('click', (ev) => { ev.stopPropagation(); loadLigand(r.smiles, r.name); });
      card.addEventListener('click', () => loadLigand(r.smiles, r.name));
      results.appendChild(card);
    });
    state.lastLigandResults = data.results;
    if (typeof simBadgeResults === 'function') simBadgeResults(data.results);

  } catch (err) {
    results.innerHTML = `<div class="browser-placeholder"><div style="color:var(--red)">Search failed: ${err.message}</div></div>`;
  }
}

function loadLigand(smiles, name) {
  if (!smiles || !smiles.trim()) {
    showToast(`${name}: PubChem returned no SMILES for this compound — try another result.`, 'error', '⚠️');
    return;
  }
  setSMILES(smiles);
  showToast(`Loaded ${name} — switch to Dock tab to begin!`, 'success', '💊');
  // Switch to dock tab
  document.querySelector('[data-tab="dock"]').click();
}


// BEGIN RELATED_PDB_UI_V1
(function () {
  'use strict';
  if (window.__relatedPdbUiV1) return;
  window.__relatedPdbUiV1 = true;

  function renderRelatedPdbs(data) {
    const grid = document.getElementById('target-grid');
    if (!grid || !grid.parentElement) return;

    let card = document.getElementById('related-pdb-card');
    if (!card) {
      card = document.createElement('div');
      card.id = 'related-pdb-card';
      card.className = 'card';
      card.style.margin = '12px 0';
      grid.parentElement.insertBefore(card, grid);
    }

    const rows = data.results || [];
    card.innerHTML = `
      <div class="card-title">🧬 Similar PDB structures</div>
      <div class="text-xs" style="color:var(--text-2);margin:5px 0 9px">
        Structures related to ${escapeHtml(data.pdb_id)}. Binding molecules are listed from RCSB.
      </div>
      ${rows.length ? rows.map(row => `
        <div style="margin:7px 0;padding:8px;border:1px solid var(--border);border-radius:8px">
          <div style="display:flex;justify-content:space-between;gap:8px;align-items:flex-start">
            <div class="text-xs">
              <b>${escapeHtml(row.pdb_id)}</b>
              ${row.resolution != null ? ` · ${escapeHtml(row.resolution)} Å` : ''}
              <br><span style="color:var(--text-2)">${escapeHtml(row.title)}</span>
              <br><span style="color:var(--text-3)">Method: ${escapeHtml(row.method || 'not listed')}</span>
              <br><span style="color:var(--text-3)">Binding molecules: ${
                (row.binding_ligands || []).length
                  ? row.binding_ligands.map(l => `${escapeHtml(l.ccd_id)} (${escapeHtml(l.name)})`).join(' · ')
                  : 'none listed'
              }</span>
            </div>
            <button class="btn btn-sm" data-related-pdb="${escapeHtml(row.pdb_id)}">Use PDB</button>
          </div>
        </div>
      `).join('') : '<div class="text-xs" style="color:var(--text-3)">No related experimental structures were returned.</div>'}`;

    card.querySelectorAll('[data-related-pdb]').forEach(button => {
      button.addEventListener('click', () => window.importPdb(button.dataset.relatedPdb));
    });
  }

  const originalImportPdb = window.importPdb;
  window.importPdb = async function (pdbId) {
    const chosen = (pdbId || (document.getElementById('pdb-import-input') || {}).value || '').trim().toUpperCase();
    const data = await originalImportPdb.apply(this, arguments);

    if (data && chosen) {
      try {
        const related = await apiFetch(`/api/pdb/${encodeURIComponent(chosen)}/related?limit=6`);
        renderRelatedPdbs(related);
      } catch (err) {
        console.warn('Related PDB search failed:', err);
      }
    }
    return data;
  };
})();
// END RELATED_PDB_UI_V1
