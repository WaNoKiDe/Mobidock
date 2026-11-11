// ══ CELL 43 · static/js/literature.js ══
/* ─── literature.js — literature — Europe PMC papers for the selected target ─── */


// ── Literature ────────────────────────────────────────────────────────────────
let literatureRequest = 0;
async function loadLiterature(targetId) {
  const card = document.getElementById('literature-card');
  const content = document.getElementById('literature-content');
  const titleName = document.getElementById('lit-target-name');
  if (!card || !content || !titleName || !targetId) return;
  const seq = ++literatureRequest;
  card.style.display = 'block';
  content.innerHTML = '<div style="text-align:center;color:var(--text-3)">Searching Europe PMC... ⏳</div>';
  try {
    const res = await fetch('/api/literature?target_id=' + targetId);
    if (!res.ok) throw new Error(`literature service returned ${res.status}`);
    const data = await res.json();
    if (seq !== literatureRequest || state.selectedTarget !== targetId) return;
    titleName.textContent = data.target || (state.targets[targetId] || {}).name || targetId;
    if (!data.papers || data.papers.length === 0) {
      content.innerHTML = '<div style="color:var(--text-3)">No recent papers found.</div>';
      return;
    }
    let html = '<ul style="padding-left:16px; margin:0">';
    data.papers.forEach(p => {
      const link = p.doi ? `https://doi.org/${p.doi}` : `https://pubmed.ncbi.nlm.nih.gov/${p.pmid}`;
      html += `<li style="margin-bottom:8px">
        <a href="${link}" target="_blank" style="color:var(--accent-2);text-decoration:none;font-weight:bold">${p.title}</a><br>
        <span style="color:var(--text-3);font-size:11px">${p.authors} (${p.year}) - ${p.journal}</span>
      </li>`;
    });
    html += '</ul>';
    content.innerHTML = html;
  } catch (err) {
    if (seq !== literatureRequest || state.selectedTarget !== targetId) return;
    titleName.textContent = (state.targets[targetId] || {}).name || targetId;
    content.innerHTML = `<div style="color:var(--text-3)">Literature search is unavailable right now${err.message ? `: ${escapeHtml(err.message)}` : '.'}</div>`;
  }
}
