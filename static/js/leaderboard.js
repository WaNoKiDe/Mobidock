// ══ CELL 42 · static/js/leaderboard.js ══
/* ─── leaderboard.js — leaderboard — tabs and table ─── */


// ── Leaderboard ───────────────────────────────────────────────────────────────
function renderLeaderboardTabs() {
  const wrap = document.getElementById('lb-tabs');
  wrap.innerHTML = '';
  Object.keys(state.targets).forEach(tid => {
    const tab = document.createElement('div');
    tab.className = `lb-tab${tid === state.lbCurrentTab ? ' active' : ''}`;
    tab.textContent = tid;
    tab.setAttribute('id', `lb-tab-${tid}`);
    tab.addEventListener('click', () => {
      state.lbCurrentTab = tid;
      document.querySelectorAll('.lb-tab').forEach(t => t.classList.remove('active'));
      tab.classList.add('active');
      refreshLeaderboard();
    });
    wrap.appendChild(tab);
  });
}

async function refreshLeaderboard() {
  try {
    const data = await apiFetch(`/api/leaderboard?target_id=${state.lbCurrentTab}&top_n=10`);
    const entries = data[state.lbCurrentTab] || [];
    const tbody = document.getElementById('lb-body');
    tbody.innerHTML = '';

    if (!entries.length) {
      tbody.innerHTML = '<tr><td colspan="6" class="lb-empty">No scores yet — be the first to dock! 🚀</td></tr>';
      return;
    }

    entries.forEach((e, i) => {
      const rank = i + 1;
      const rankClass = rank === 1 ? 'gold' : rank === 2 ? 'silver' : rank === 3 ? 'bronze' : '';
      const rankIcon = rank === 1 ? '🥇' : rank === 2 ? '🥈' : rank === 3 ? '🥉' : `#${rank}`;
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td class="lb-rank ${rankClass}">${rankIcon}</td>
        <td>${escapeHtml(e.player)}</td>
        <td class="lb-energy">${e.energy.toFixed(2)}</td>
        <td class="lb-qed">${e.qed.toFixed(3)}</td>
        <td class="lb-synth" title="SA score: 1 = easy to make, 10 = very hard">${(e.synth === undefined || e.synth === null) ? '—' : Number(e.synth).toFixed(1)}</td>
        <td style="color:var(--accent);font-weight:600;font-family:'JetBrains Mono',monospace">${e.score.toFixed(2)}</td>
      `;
      tbody.appendChild(tr);
    });
  } catch (e) {
    console.warn('Leaderboard fetch failed:', e);
  }
}
