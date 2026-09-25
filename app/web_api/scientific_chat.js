"use strict";
const base = '/api/v1/scientific';
const $ = id => document.getElementById(id);
const runningStates = new Set(['queued', 'preparing', 'running']);
let token = '', identity = null, active = null, conversation = null, displayed = '';
let loaded = false, submitting = false, polling = false, pollAgain = false;
let pollTimer = null, navigation = 0, stopping = null, sidebarOpen = false, activityVersion = 0;
const drafts = new Map();
function element(tag, text = '', className = '') {
  const node = document.createElement(tag);
  node.textContent = text;
  if (className) node.className = className;
  return node;
}
function formattedMessage(text,presentation){
  const body=element('div','','answer');
  // Only server-generated Markdown HTML enters this sink. Raw model HTML is disabled.
  if(presentation?.html!==undefined){body.innerHTML=presentation.html;body.querySelectorAll('table').forEach(table=>{const scroll=element('div','','table-scroll');scroll.tabIndex=0;scroll.setAttribute('role','region');scroll.setAttribute('aria-label','Scrollable answer table');table.replaceWith(scroll);scroll.append(table)})}
  else{body.textContent=text;body.style.whiteSpace='pre-wrap'}
  return body;
}
function appendStructures(container,presentation){
  if(!presentation?.structures?.length)return;
  container.append(element('div','Structures from the SMILES in this message','structure-heading'));
  const gallery=element('div','','structure-gallery');
  presentation.structures.forEach((structure,index)=>{
    const figure=element('figure','','structure');const image=element('img');image.src=structure.image_url;image.alt='Molecular structure: '+structure.smiles;image.loading='lazy';
    const caption=element('figcaption','Structure '+(index+1));caption.append(element('code',structure.smiles));
    const download=element('a','Download SVG');download.href=structure.image_url;download.download='structure-'+(index+1)+'.svg';caption.append(download);figure.append(image,caption);gallery.append(figure);
  });container.append(gallery);
}
const basisLabels={reported:'Reported · agent attributed',computed:'Computed · recorded execution',proposed:'Proposed · hypothesis',unknown:'Unknown',input:'User input'};
function showStructured(view){
  const section=element('section','','scientific-view');
  if(view.error){section.append(element('p',view.error,'scientific-note'));return section}
  const sources=new Map(view.sources.map(source=>[source.id,source]));const molecules=new Map(view.molecules.map(molecule=>[molecule.id,molecule]));const steps=new Map(view.steps.map(step=>[step.id,step]));
  function links(container,ids){for(const id of ids){const source=sources.get(id);if(!source)continue;const link=element('a',source.title);link.href=source.artifact_url;link.target='_blank';link.rel='noopener';link.title=source.locator;container.append(link,document.createTextNode(' '))}}
  function attributed(container,item){container.append(element('span',basisLabels[item.basis]||item.basis,'basis basis-'+item.basis));const citations=element('span','','evidence');links(citations,item.source_ids);container.append(citations);for(const limitation of item.limitations)container.append(element('p',limitation,'scientific-note'))}
  function claim(item){const box=element('div','','attributed-claim');box.append(element('div',item.text));attributed(box,item);return box}
  if(view.molecules.length||view.steps.length||view.claims.length)section.append(element('p','Structured scientific view · agent authored, not independently reviewed. Drawings depict declared structures; they do not validate feasibility.','scientific-note'));
  if(view.molecules.length){section.append(element('h3','Molecules and intermediates'));const gallery=element('div','','structure-gallery');for(const molecule of view.molecules){const card=element('figure','','structure');card.append(element('figcaption',molecule.name+(view.target_molecule_ids.includes(molecule.id)?' · target':'')));if(molecule.image_url){const image=element('img');image.src=molecule.image_url;image.alt=molecule.name+' · '+molecule.smiles;image.loading='lazy';card.append(image);const download=element('a','Download SVG');download.href=molecule.image_url;download.download=molecule.id+'.svg';card.append(download)}else card.append(element('p','Structure could not be drawn. Original notation is retained.','scientific-note'));card.append(element('code',molecule.smiles));attributed(card,molecule);gallery.append(card)}section.append(gallery)}
  function stepCard(step){const card=element('article','','scientific-step');card.append(element('h4',step.id+' · '+step.title));attributed(card,step);card.append(element('p',step.reactant_ids.map(id=>molecules.get(id).name).join(' + ')+' → '+step.product_ids.map(id=>molecules.get(id).name).join(' + ')));
    if(step.image_url){const image=element('img','','reaction-scheme');image.src=step.image_url;image.alt='Declared reaction scheme for '+step.title;image.loading='lazy';card.append(image);const download=element('a','Download reaction SVG');download.href=step.image_url;download.download=step.id+'-reaction.svg';card.append(download)}else card.append(element('p','Reaction scheme could not be drawn; inspect the declared SMILES.','scientific-note'));
    const notation=element('details');notation.append(element('summary','Reaction SMILES'),element('code',step.reaction_smiles));card.append(notation);
    if(step.after_step_ids.length)card.append(element('p','Depends on: '+step.after_step_ids.join(', '),'muted'));
    card.append(element('h4','Conditions'));if(step.conditions.length)step.conditions.forEach(item=>card.append(claim(item)));else card.append(element('p','No condition details supplied.','muted'));
    card.append(element('h4','Yield'));card.append(step.yield_info?claim(step.yield_info):element('p','Not supplied.','muted'));return card;
  }
  const routed=new Set();for(const route of view.routes){const box=element('section');box.append(element('h3',route.title));const diagram=element('img','','route-overview');diagram.src=route.image_url;diagram.alt='Step dependencies for '+route.title;box.append(diagram);for(const limitation of route.limitations)box.append(element('p',limitation,'scientific-note'));if(route.unreached_target_ids.length)box.append(element('p','Incomplete route: no terminal product with the declared target ID(s) '+route.unreached_target_ids.join(', ')+'.','scientific-note'));for(const id of route.step_ids){box.append(stepCard(steps.get(id)));routed.add(id)}section.append(box)}
  const standalone=view.steps.filter(step=>!routed.has(step.id));if(standalone.length){section.append(element('h3','Reaction steps'));standalone.forEach(step=>section.append(stepCard(step)))}
  if(view.claims.length){section.append(element('h3','Findings and proposals'));view.claims.forEach(item=>section.append(claim(item)))}
  return section;
}

async function request(path, options = {}) {
  const response = await fetch(base + path, {
    ...options, headers: {'Content-Type': 'application/json', 'X-Scientific-Token': token, ...options.headers},
  });
  const data = await response.json();
  if (!response.ok) throw new Error(typeof data.detail === 'string' ? data.detail : JSON.stringify(data.detail));
  return data;
}

function disclosure(title, key, content) {
  const box = element('details', '', 'disclosure');
  box.dataset.key = key;
  box.append(element('summary', title), content);
  return box;
}

function sourcesPanel(turn) {
  const container = element('div');
  const sources = turn.structured_presentation?.sources || [];
  const linked = new Set();
  function link(row, title, href) {
    const a = element('a', title);
    a.href = href; a.target = '_blank'; a.rel = 'noopener noreferrer'; row.append(a);
  }
  for (const source of sources) {
    const row = element('div', '', 'source-card');
    row.append(element('strong', source.title), element('p', source.locator, 'muted'));
    link(row, 'Inspect saved evidence', source.artifact_url);
    if (source.url) link(row, 'Original source ↗', source.url);
    container.append(row); linked.add(source.artifact_ref);
  }
  for (const [index, ref] of (turn.answer.evidence_refs || []).entries()) {
    if (linked.has(ref)) continue;
    const row = element('div', '', 'source-card');
    link(row, 'Recorded evidence ' + (index + 1), base + '/conversations/' + identity + '/artifacts/' + encodeURIComponent(ref));
    container.append(row);
  }
  return container;
}

function answerCard(turn) {
  const card = element('article', '', 'assistant');
  if (turn.answer) {
    const answer = turn.answer;
    card.append(formattedMessage(answer.answer_markdown, turn.answer_presentation));
    const view = turn.structured_presentation;
    if (view && (view.error || view.molecules?.length || view.steps?.length || view.claims?.length)) {
      const label = view.error ? 'Scientific details unavailable' : 'Structures & scientific details' +
        (view.molecules.length ? ' · ' + view.molecules.length + ' molecules' : '');
      card.append(disclosure(label, turn.id + ':structures', showStructured(view)));
    } else if (turn.answer_presentation?.structures?.length) {
      const structures = element('div'); appendStructures(structures, turn.answer_presentation);
      card.append(disclosure('View structures', turn.id + ':structures', structures));
    }
    if (answer.evidence_refs?.length || view?.sources?.length) {
      card.append(disclosure('Sources & evidence', turn.id + ':sources', sourcesPanel(turn)));
    }
    if (answer.uncertainties?.length) {
      const notes = element('div'); answer.uncertainties.forEach(text => notes.append(element('p', text)));
      card.append(disclosure('Uncertainty & missing information', turn.id + ':uncertainty', notes));
    }
    if (answer.needs_user_input) card.append(element('p', 'Add the requested information below to continue.', 'muted'));
    const actions = element('div', '', 'answer-actions');
    const copy = element('button', 'Copy', 'copy-button'); copy.type = 'button';
    copy.setAttribute('aria-label', 'Copy answer');
    copy.onclick = async () => {
      try { await navigator.clipboard.writeText(answer.answer_markdown); copy.textContent = 'Copied'; }
      catch { copy.textContent = 'Unable to copy'; }
      setTimeout(() => { copy.textContent = 'Copy'; }, 2000);
    };
    actions.append(copy, element('span', 'Research answer · not independently reviewed'));
    card.append(actions);
  } else if (runningStates.has(turn.status)) {
    const progress = element('div', '', 'progress'); progress.id = 'live-progress';
    const line = element('div', '', 'progress-line');
    const label = element('span', 'Starting investigation…'); label.id = 'progress-label'; label.setAttribute('role', 'status');
    const elapsed = element('span', '', 'elapsed'); elapsed.id = 'elapsed';
    const spinner = element('span', '', 'spinner'); spinner.setAttribute('aria-hidden', 'true');
    line.append(spinner, label, elapsed);
    const detail = element('details'); detail.dataset.key = turn.id + ':activity';
    detail.append(element('summary', 'View activity'));
    const events = element('ol', '', 'activity-list'); events.id = 'activity-list'; detail.append(events);
    progress.append(line, element('p', 'The answer will appear here when it is ready.', 'progress-caption'), detail);
    card.append(progress);
  } else {
    const stopped = turn.status === 'cancelled';
    card.append(element('p', stopped ? 'Stopped. Your question and any saved evidence are preserved.' :
      'This investigation could not finish.', stopped ? 'muted' : 'error-message'));
    if (turn.error && !stopped) card.append(disclosure('Show details', turn.id + ':error', element('p', turn.error.message)));
  }
  return card;
}

function nearBottom() {
  const scroller = $('scroll-area');
  return scroller.scrollHeight - scroller.scrollTop - scroller.clientHeight < 160;
}
function toBottom() { $('scroll-area').scrollTop = $('scroll-area').scrollHeight; }
function renderConversation(data) {
  const firstLoad = conversation === null;
  conversation = data;
  document.title = data.title + ' · Scientific workspace';
  const signature = JSON.stringify(data.turns.map(t => [t.id, t.status, t.answer_ref, t.error]));
  if (signature !== displayed) {
    const follow = firstLoad || nearBottom();
    const open = new Set(Array.from($('messages').querySelectorAll('details[open]')).map(node => node.dataset.key).filter(Boolean));
    $('messages').replaceChildren();
    for (const turn of data.turns) {
      const question = element('div', '', 'user');
      const bubble = element('div', '', 'user-bubble'); bubble.append(formattedMessage(turn.question, turn.question_presentation));
      question.append(bubble);
      if (turn.question_presentation?.structures?.length) {
        const structures = element('div'); appendStructures(structures,turn.question_presentation);
        question.append(disclosure('View structure', turn.id + ':input', structures));
      }
      $('messages').append(question, answerCard(turn));
    }
    $('messages').querySelectorAll('details').forEach(node => { if (open.has(node.dataset.key)) node.open = true; });
    displayed = signature;
    if (follow) requestAnimationFrame(toBottom);
  }
  updateProgress();
}

function elapsedText(start) {
  const timestamp = Date.parse(start);
  if (!Number.isFinite(timestamp)) return '';
  const seconds = Math.max(0, Math.floor((Date.now() - timestamp) / 1000));
  return seconds >= 60 ? Math.floor(seconds / 60) + 'm ' + seconds % 60 + 's' : seconds + 's';
}
function activityLabel(event) {
  const names = {command_execution: 'Analysis tool', web_search: 'Web search', mcp_tool_call: 'Scientific tool',
    file_change: 'Saving investigation files', todo_list: 'Updating the plan'};
  const label = names[event.kind] || 'Agent activity';
  return label + (['completed', 'item.completed'].includes(event.status) ? ' · finished' :
    ['failed', 'item.failed'].includes(event.status) ? ' · failed' : ' · running');
}
function updateProgress() {
  const turn = active?.conversation_id === identity ? active : conversation?.turns.at(-1);
  const label = $('progress-label');
  if (label && turn && runningStates.has(turn.status)) {
    label.textContent = stopping === identity ? 'Stopping investigation…' :
      turn.status === 'preparing' ? 'Preparing your workspace…' : turn.status === 'queued' ? 'Starting investigation…' :
      turn.repair_attempts ? 'Refining the answer…' : 'Investigating…';
    $('elapsed').textContent = elapsedText(turn.created_at);
    const events = $('activity-list');
    const signature = JSON.stringify(turn.progress || []);
    if (events.dataset.signature !== signature) {
      events.replaceChildren();
      for (const event of (turn.progress || []).slice(-8)) {
        const row = element('li'); const at = new Date(event.at);
        row.append(element('time', Number.isNaN(at.getTime()) ? '' : at.toLocaleTimeString([], {hour:'2-digit', minute:'2-digit', second:'2-digit'})),
          element('span', activityLabel(event))); events.append(row);
      }
      if (!events.children.length) events.append(element('li', 'Waiting for the first activity update…'));
      events.dataset.signature = signature;
    }
  }
  if (active) $('active-text').textContent = (stopping ? 'Stopping investigation' : 'An investigation is running') + ' · ' + elapsedText(active.created_at);
}
function resizeInput() {
  $('question').style.height = 'auto'; $('question').style.height = Math.min($('question').scrollHeight, 180) + 'px';
}
function updateControls() {
  $('send').disabled = !loaded || submitting || !!active || !$('question').value.trim();
  $('send').hidden = !!active; $('cancel').hidden = !active;
  $('cancel').disabled = !!stopping; $('cancel').setAttribute('aria-label', stopping ? 'Stopping investigation' : 'Stop investigation');
  $('new').disabled = submitting; $('question').disabled = submitting;
  $('global-activity').hidden = !active || active.conversation_id === identity;
  $('stop-active').disabled = !!stopping;
  $('composer-hint').textContent = submitting ? 'Sending…' : !loaded ? 'Connecting…' :
    active ? (stopping ? 'Stopping…' : 'You can write a follow-up while the agent works') : 'Enter to send · Shift + Enter for a new line';
  updateProgress();
}

async function list() {
  const rows = await request('/conversations'); $('saved').replaceChildren();
  for (const row of rows) {
    const button = element('button', '', 'saved'); button.title = row.title;
    button.setAttribute('aria-current', identity === row.id ? 'page' : 'false');
    if (active?.conversation_id === row.id) { const dot = element('span', '', 'running-dot'); dot.setAttribute('aria-label', 'Running'); button.append(dot); }
    button.append(element('span', row.title, 'saved-title'));
    button.onclick = () => { if (!submitting) openConversation(row.id); }; $('saved').append(button);
  }
  if (!rows.length) $('saved').append(element('p', 'Your chats will appear here.', 'nav-label'));
}
async function loadConversation(id, epoch) {
  const data = await request('/conversations/' + id);
  if (epoch === navigation && id === identity) renderConversation(data);
}
async function poll() {
  clearTimeout(pollTimer);
  if (polling) { pollAgain = true; return; }
  polling = true;
  try {
    const previous = active;
    const version = activityVersion;
    const activity = await request('/activity');
    // A response begun before a submitted turn must not hide its Stop button.
    if (version === activityVersion) active = activity.active;
    if (stopping && active?.conversation_id !== stopping) stopping = null;
    updateControls();
    if (identity && (!conversation || active?.conversation_id === identity || previous?.conversation_id === identity ||
        runningStates.has(conversation.turns.at(-1)?.status))) await loadConversation(identity, navigation);
    if (previous?.conversation_id !== active?.conversation_id) await list();
    if ($('error').dataset.network === 'true') { $('error').textContent = ''; delete $('error').dataset.network; }
  } catch (error) {
    $('error').textContent = 'Connection interrupted. Retrying… ' + error.message;
    $('error').dataset.network = 'true';
  } finally {
    polling = false;
    pollTimer = setTimeout(poll, pollAgain ? 0 : active ? 1500 : 4000); pollAgain = false;
  }
}

function rememberDraft() { drafts.set(identity || 'new', $('question').value); }
function setSidebar(open) {
  sidebarOpen = open; $('app').classList.toggle('sidebar-closed', !open);
  $('sidebar').inert = !open; $('menu').setAttribute('aria-expanded', String(open));
}
async function openConversation(id) {
  rememberDraft(); identity = id; navigation++; conversation = null; displayed = '';
  const epoch = navigation; location.hash = id; $('intro').hidden = true;
  $('messages').replaceChildren(element('p', 'Loading conversation…', 'muted'));
  $('error').textContent = ''; $('question').value = drafts.get(id) || ''; resizeInput(); updateControls();
  if (matchMedia('(max-width:760px)').matches) setSidebar(false);
  try { await loadConversation(id, epoch); await list(); await poll(); }
  catch (error) { if (epoch === navigation) $('error').textContent = error.message; }
}
function newChat() {
  if (submitting) return;
  rememberDraft(); identity = null; navigation++; conversation = null; displayed = '';
  location.hash = ''; document.title = 'Scientific workspace'; $('messages').replaceChildren();
  $('intro').hidden = false; $('error').textContent = ''; $('jump-bottom').hidden = true;
  $('question').value = drafts.get('new') || ''; resizeInput(); updateControls();
  if (matchMedia('(max-width:760px)').matches) setSidebar(false);
  list().catch(error => { $('error').textContent = error.message; }); $('question').focus();
}
async function stopActive() {
  const target = active?.conversation_id;
  if (!target || stopping) return;
  stopping = target; updateControls(); $('error').textContent = '';
  try {
    const result = await request('/conversations/' + target + '/cancel', {method:'POST'});
    if (!result.cancellation_requested) stopping = null;
    await poll();
  } catch (error) { stopping = null; $('error').textContent = error.message; }
  updateControls();
}

$('new').onclick = newChat;
$('menu').onclick = () => setSidebar(!sidebarOpen);
$('backdrop').onclick = () => { setSidebar(false); $('menu').focus(); };
document.addEventListener('keydown', event => { if (event.key === 'Escape' && matchMedia('(max-width:760px)').matches) { setSidebar(false); $('menu').focus(); } });
$('open-active').onclick = () => { if (active) openConversation(active.conversation_id); };
$('cancel').onclick = stopActive; $('stop-active').onclick = stopActive;
$('question').oninput = () => { resizeInput(); updateControls(); };
$('question').onkeydown = event => {
  if (event.key === 'Enter' && !event.shiftKey && !event.isComposing) {
    event.preventDefault(); if (!$('send').disabled) $('form').requestSubmit();
  }
};
document.querySelectorAll('[data-example]').forEach(button => { button.onclick = () => { $('question').value = button.dataset.example; resizeInput(); updateControls(); $('question').focus(); }; });
$('scroll-area').onscroll = () => { $('jump-bottom').hidden = !identity || nearBottom(); };
$('jump-bottom').onclick = toBottom;
window.addEventListener('hashchange', () => {
  const id = location.hash.slice(1);
  if (id === (identity || '')) return;
  if (/^[0-9a-f]{32}$/.test(id)) openConversation(id); else newChat();
});
$('form').onsubmit = async event => {
  event.preventDefault(); const question = $('question').value.trim();
  if (!question || !loaded || active || submitting) return;
  const source = identity; submitting = true; activityVersion++; updateControls(); $('error').textContent = '';
  try {
    const turn = await request('/turns', {method:'POST', body:JSON.stringify({question, conversation_id:source})});
    $('question').value = ''; drafts.delete(source || 'new');
    activityVersion++;
    active = {conversation_id:turn.conversation_id, id:turn.turn_id, status:'queued', created_at:new Date().toISOString(), progress:[]};
    submitting = false; await openConversation(turn.conversation_id);
  } catch (error) {
    $('error').textContent = error.message; submitting = false;
    await poll(); // Discover and expose the running chat after a cross-tab race.
  }
  updateControls(); $('question').focus();
};

async function initialize() {
  setSidebar(!matchMedia('(max-width:760px)').matches);
  try {
    const config = await request('/config'); token = config.token;
    $('environment').replaceChildren(element('p', 'Agent: ' + config.runtime), element('p', 'Model: ' + config.model),
      element('p', 'Questions and tool context may be sent to your configured model provider. Investigation records stay local.'));
    active = (await request('/activity')).active; await list(); loaded = true; updateControls();
    const id = location.hash.slice(1);
    if (/^[0-9a-f]{32}$/.test(id)) await openConversation(id);
    else if (active) await openConversation(active.conversation_id);
    else await poll();
  } catch (error) { $('error').textContent = 'Unable to connect. ' + error.message; updateControls(); }
}
setInterval(updateProgress, 1000);
initialize();
