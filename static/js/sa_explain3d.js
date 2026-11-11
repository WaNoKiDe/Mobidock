/* Shared on-demand SA explanation. Uses one viewer to avoid WebGL context buildup. */
(function () {
  'use strict';
  if(window.__SA3D_INLINE_V2)return;window.__SA3D_INLINE_V2=true;
  const esc = x => String(x == null ? '' : x).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const n = x => Number(x).toFixed(4);
  let viewer=null, dialog=null, request=null, generation=0, observer=null;
  function source(s) { return s && (s.sa_smiles || s.smiles) || ''; }
  function button(s) {
    const smi=source(s);
    return `<div class="sa3d-inline" data-sa-smiles="${esc(smi)}" style="margin-top:12px;min-height:100px">${smi?'Loading 3D ligand and Morgan contributions…':'Analyse the molecule again after restarting the server to load its 3D explanation.'}</div>`;
  }
  window.synthDriversHTML = button;
  window.plSynthHTML = function(s) {
    let h='<div class="gd-sec pl-tile"><div class="gd-h">🧪 Synthesizability <span>structural SA estimate</span></div>';
    if(!s||!s.available)return h+'<p>SA scorer unavailable. Check the server setup.</p></div>';
    h+=`<div style="display:flex;align-items:baseline;gap:8px;flex-wrap:wrap"><span class="pl-big" style="color:${esc(s.colour)}">${esc(s.score)}</span><span>/100</span><span class="gd-chip">${esc(s.tier_label)} · SA ${esc(s.sa)}</span></div>`;
    if(typeof window.synthMeterHTML==='function')h+=window.synthMeterHTML(s.sa);
    return h+button(s)+'</div>';
  };
  function close() {
    generation++;if(request)request.abort();request=null;
    if(observer)observer.disconnect();observer=null;
    if(viewer){viewer.clear();viewer=null;}
    if(dialog){dialog.remove();dialog=null;}
  }
  function makeDialog() {
    close();dialog=document.createElement('dialog');
    dialog.style.cssText='width:min(1100px,95vw);max-height:92vh;overflow:auto;background:#101827;color:#e5e7eb;border:1px solid #475569;border-radius:14px;padding:18px;';
    dialog.innerHTML='<button type="button" data-sa-close style="float:right">Close ×</button><h2>SA score · Morgan environments on your ligand</h2><div class="sa-body">Loading fresh fragment data…</div>';
    dialog.addEventListener('cancel',e=>{e.preventDefault();close();});
    dialog.querySelector('[data-sa-close]').onclick=close;
    document.body.appendChild(dialog);dialog.showModal();return dialog.querySelector('.sa-body');
  }
  async function open(smi) {
    const body=makeDialog(),token=generation;
    request=new AbortController();const timer=setTimeout(()=>request&&request.abort(),60000);
    try{
      const response=await fetch('/api/synth/explain3d',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({smiles:smi,mode:'off'}),signal:request.signal,cache:'no-store'});
      if(!response.ok){let detail='';try{detail=(await response.json()).detail;}catch(_){}throw new Error(detail||`Server returned ${response.status}. Run the corrective cell and restart the server.`);}
      const e=await response.json();if(token!==generation||!dialog)return;
      if(!window.$3Dmol&&typeof window.loadViewerLibrary==='function')await window.loadViewerLibrary();
      if(e.version!=='sa3d-1')throw new Error('Unexpected explanation version. Restart the patched server.');
      render(body,e);
    }catch(err){if(token===generation&&dialog)body.textContent=err.name==='AbortError'?'Request timed out. Reopen to retry; inspect the server log if this repeats.':err.message;}
    finally{clearTimeout(timer);}
  }
  function render(body,e) {
    let viewer=null, observer=null;
    const fragmentName=r=>[...new Set((r.environments||[]).map(x=>{const atom=e.elements&&e.elements[x.centre]||'atom';const desc=x.fragment_label||`${atom} neighborhood · radius ${x.radius}`;return x.fragment_smiles?`${desc} · ${x.fragment_smiles}`:desc;}))].join(' / ')||'unlabeled Morgan environment';
    const rows=e.rows.map((r,i)=>`<tr><td><button type="button" data-sa-row="${i}" title="Morgan hash ${esc(r.id)}">${esc(fragmentName(r))}<br><small>Morgan ID ${esc(r.id)}</small></button></td><td>${r.count}</td><td>${r.matched?'Matched':'Unseen'}</td><td>${n(r.weight)}</td><td>${n(r.contribution)}</td></tr>`).join('');
    body.innerHTML=`<p><b>SA ${e.sa} · ${e.score}/100</b> · ${e.matched_unique}/${e.unique} unique Morgan IDs matched (${e.matched_occurrences}/${e.occurrences} occurrences).</p>
      <p style="overflow-wrap:anywhere">Ligand: ${esc(e.smiles)}</p>
      <div style="display:flex;flex-wrap:wrap;gap:16px"><section style="flex:1 1 100%;min-width:0"><label>Local fragment structure (SMILES) · Morgan ID <select data-sa-id style="width:100%">${e.rows.map((r,i)=>`<option value="${i}">${fragmentName(r)} · Morgan ID ${r.id} · ${r.matched?'matched':'unseen'} · count ${r.count}</option>`).join('')}</select></label><label>Occurrence <select data-sa-occ style="width:100%"></select></label>
      <div data-sa-view style="height:300px;width:100%;position:relative;background:#07111f;border-radius:10px;margin-top:8px"></div>
      <p data-sa-selection></p><button type="button" data-sa-reset>Reset view</button><p>${esc(e.viewer_note)}</p><p>Cyan: selected fragment atoms. Gold: centre atom(s). Grey: remaining ligand. Drag to rotate; scroll to zoom.</p></section>
      <section style="flex:1 1 100%;min-width:0"><p>Each fragment is shown as a local SMILES structure; its Morgan hash is listed alongside for technical lookup. Click a fragment to highlight it. Select an occurrence, or all occurrences of the same ID. Atom indices are zero-based.</p><div style="overflow:auto;max-height:430px"><table class="gd-table"><thead><tr><th>Local fragment (SMILES) · Morgan ID</th><th>Count</th><th>Lookup</th><th>Weight</th><th>Weighted term</th></tr></thead><tbody>${rows}</tbody></table></div></section></div>
      <h3>How the score is calculated</h3><p>Fragment term = Σ(count × table weight) / ${e.occurrences} = <b>${n(e.fragment_mean)}</b>. An unseen ID uses weight −4. Weighted terms are in raw-score units, not final SA points.</p>
      <p>Complexity penalties: ${Object.entries(e.penalties).map(([k,v])=>`${esc(k)} ${n(v)}`).join(' · ')}.<br>Density/symmetry correction: +${n(e.density_correction)}.</p>
      <p>Raw = ${n(e.fragment_mean)} − ${n(e.penalty_total)} + ${n(e.density_correction)} = <b>${n(e.raw)}</b>.<br>Mapped SA = 11 − (raw + 5) × 9 / 6.5 = ${n(e.mapped)}. For mapped SA &gt; 8 apply 8 + ln(mapped SA − 8), then clamp to 1–10. Final SA: <b>${e.sa}</b>.</p>
      <p>Game display: normalize the rounded SA with (10 − SA)/9, round to four decimals, then multiply by 100 and round: ${e.score}/100. This is not a probability of successful synthesis.</p>
      <p>Verified against the installed RDKit scorer. Morgan environments overlap (radii 0–2 bonds). A match means the hashed ID exists in the contribution table; it is not a reference-drug similarity match or proof of a synthesis route.</p>`;
    const box=body.querySelector('[data-sa-view]'),pick=body.querySelector('[data-sa-id]'),occ=body.querySelector('[data-sa-occ]'),label=body.querySelector('[data-sa-selection]');
    if(e.molblock&&window.$3Dmol){
      try{
        viewer=window.$3Dmol.createViewer(box,{backgroundColor:'#07111f',antialias:true});
        const model=viewer.addModel(e.molblock,'sdf');
        const atoms=model.selectedAtoms({});
        if(atoms.length!==e.atom_count||atoms.some((a,i)=>a.elem!==e.elements[i]))throw new Error('Viewer atom mapping mismatch; highlights disabled.');
        viewer.zoomTo();viewer.render();
        observer=new ResizeObserver(()=>{if(viewer){viewer.resize();viewer.render();}});observer.observe(box);
      }catch(err){if(viewer)viewer.clear();viewer=null;box.textContent=err.message;}
    }else box.textContent=e.molblock?'3Dmol.js did not load. Check your internet connection and reload the page. The calculation below is still available.':e.viewer_note;
    function highlight(){
      const row=e.rows[Number(pick.value)];
      const envs=occ.value==='all'?row.environments:[row.environments[Number(occ.value)]];
      const ids=[...new Set(envs.flatMap(x=>x.atoms))],centres=[...new Set(envs.map(x=>x.centre))];
      label.textContent=`Local fragment ${fragmentName(row)} · Morgan ID ${row.id} · ${row.matched?'matched table entry':'unseen: default weight −4'} · atoms ${ids.join(', ')} · weight ${n(row.weight)} · total contribution ${n(row.contribution)}.`;
      if(viewer){viewer.setStyle({},{stick:{color:'#64748b',radius:.12},sphere:{color:'#64748b',scale:.22}});viewer.setStyle({index:ids},{stick:{color:'#22d3ee',radius:.2},sphere:{color:'#22d3ee',scale:.3}});viewer.setStyle({index:centres},{stick:{color:'#fbbf24',radius:.2},sphere:{color:'#fbbf24',scale:.36}});viewer.render();}
    }
    function choose(){const row=e.rows[Number(pick.value)];occ.innerHTML='<option value="all">All occurrences of this ID</option>'+row.environments.map((x,i)=>`<option value="${i}">Occurrence ${i+1}: centre ${x.centre}, radius ${x.radius}</option>`).join('');occ.value='0';highlight();}
    pick.onchange=choose;occ.onchange=highlight;
    body.querySelectorAll('[data-sa-row]').forEach(b=>b.onclick=()=>{pick.value=b.dataset.saRow;choose();});
    body.querySelector('[data-sa-reset]').onclick=()=>{if(viewer){viewer.zoomTo();viewer.render();}};choose();
  }
  document.addEventListener('click',event=>{const b=event.target.closest('[data-sa3d]');if(b&&!b.disabled)open(b.dataset.sa3d);});

  window.SA3DPanel = button;
  window.SA3DScreen = window.plSynthHTML;
  const cache=new Map();
  async function hydrate(el){
    if(el.dataset.saLoaded||!el.dataset.saSmiles)return;
    el.dataset.saLoaded='1';
    const smiles=el.dataset.saSmiles;
    try{
      let pending=cache.get(smiles);
      if(!pending){
        pending=fetch('/api/synth/explain3d',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({smiles,mode:'off'}),cache:'no-store'}).then(async r=>{const e=await r.json();if(!r.ok)throw Error(e.detail||'Explanation request failed');return e;});
        cache.set(smiles,pending);if(cache.size>16)cache.delete(cache.keys().next().value);
      }
      const e=await pending;
      if(!window.$3Dmol&&typeof window.loadViewerLibrary==='function')await window.loadViewerLibrary();
      if(el.isConnected)render(el,e);
    }catch(err){cache.delete(smiles);el.textContent='3D explanation could not load: '+err.message;}
  }
  let queued=false;
  function scan(){if(queued)return;queued=true;requestAnimationFrame(()=>{queued=false;document.querySelectorAll('.sa3d-inline:not([data-sa-loaded])').forEach(el=>{if(el.getBoundingClientRect().width>0)hydrate(el);});});}
  function start(){new MutationObserver(scan).observe(document.body,{childList:true,subtree:true});scan();document.addEventListener('click',scan);}
  if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',start);else start();
})();
