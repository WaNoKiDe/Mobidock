/* Static molecular comparison cards with explicit shared and differing groups. */
function simInteractiveEscape(value) {
  return escapeHtml(String(value == null ? '' : value));
}

function simFeatureList(items, emptyText, tone) {
  if (!items || !items.length) return `<p class="sim-feature-empty">${simInteractiveEscape(emptyText)}</p>`;
  return `<ul class="sim-feature-list ${tone}">${items.map(item =>
    `<li><span>${simInteractiveEscape(item.label)}</span>${item.count > 1 ? `<small>× ${item.count}</small>` : ''}</li>`
  ).join('')}</ul>`;
}

function simFeatureSummary(row) {
  const groups = row.group_features || {};
  const shared = groups.shared || [];
  const onlyQuery = groups.only_query || [];
  const onlyCandidate = groups.only_candidate || [];
  return `<section class="sim-feature-summary" aria-label="Molecular feature comparison">
    <div class="sim-feature-summary-title">Chemical feature comparison</div>
    <div class="sim-feature-grid">
      <section><h4>Shared features</h4>${simFeatureList(shared, 'No recognized functional groups or rings are shared.', 'shared')}</section>
      <section><h4>Only in loaded ligand</h4>${simFeatureList(onlyQuery, 'No differing groups or substituents.', 'different')}</section>
      <section><h4>Only in compared molecule</h4>${simFeatureList(onlyCandidate, 'No differing groups or substituents.', 'different')}</section>
    </div>
  </section>`;
}

window.simRender = function (data, box) {
  if (!data.results || !data.results.length) {
    box.innerHTML = '<div class="text-xs" style="color:var(--text-3)">No similar ChEMBL molecules found.</div>';
    return;
  }
  const sourceLink = row => row.source === 'chembl' && row.chembl
    ? `<a href="https://www.ebi.ac.uk/chembl/compound_report_card/${encodeURIComponent(row.chembl)}/" target="_blank" rel="noopener" class="gd-chip" style="border-color:#4cc9f0;color:#4cc9f0;text-decoration:none">${simInteractiveEscape(row.chembl)}</a>`
    : '';
  const metric = data.ranked_by === 'morgan' ? 'ECFP4 Morgan' : 'MACCS-key';
  box.innerHTML = data.results.map(row => `<article class="gd-sec sim-comparison-card" style="margin-bottom:12px">
    <div class="sim-comparison-header">
      <div><strong>${simInteractiveEscape(row.name)}</strong> ${sourceLink(row)}</div>
      <div class="sim-comparison-actions"><span><b>${simPercentDisplay(row)}</b> ${metric} similarity</span><button class="btn btn-sm sim-load" data-i="${row.index}">Load ligand</button></div>
    </div>
    <div class="sim-comparison-bar" role="meter" aria-label="${metric} fingerprint Tanimoto similarity" aria-valuemin="0" aria-valuemax="100" aria-valuenow="${simPercentWidth(row)}"><div style="width:${simPercentWidth(row)}%"></div></div>
    <p class="text-xs" style="color:var(--text-3);margin:0 0 8px">Fingerprint overlap measures matching local atom patterns, not full molecular identity. ${row.fingerprint_identical && !row.same ? 'The fingerprints match, but the molecular structures differ; red highlights show the unmatched parts.' : 'Red highlights show unmatched atoms and groups.'}</p>
    <div class="sim-pair-structures">
      <figure><figcaption>Loaded ligand</figcaption>${row.query_svg || ''}</figure>
      <figure><figcaption>${simInteractiveEscape(row.name)}</figcaption>${row.candidate_svg || ''}</figure>
    </div>
    <p class="text-xs" style="color:var(--text-3);margin:0 0 8px">Green: shared structural groups · Red: groups or atoms that differ.</p>
    ${simFeatureSummary(row)}
  </article>`).join('');
  box.querySelectorAll('.sim-load').forEach(button => button.addEventListener('click', () => {
    const row = data.results.find(result => String(result.index) === button.dataset.i);
    if (row) loadLigand(row.smiles, row.name);
  }));
};
