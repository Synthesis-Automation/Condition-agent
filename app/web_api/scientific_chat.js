"use strict";
const base = '/api/v1/scientific';
const $ = id => document.getElementById(id);
const runningStates = new Set(['queued', 'preparing', 'running']);
let token = '', identity = null, active = null, conversation = null, displayed = '';
let loaded = false, submitting = false, polling = false, pollAgain = false;
let pollTimer = null, navigation = 0, stopping = null, sidebarOpen = false, activityVersion = 0;
const drafts = new Map();
function updateThemeButton() {
  const dark = document.documentElement.dataset.theme === 'dark';
  $('theme').textContent = dark ? 'Light mode' : 'Dark mode';
  $('theme').title = dark ? 'Switch to light theme' : 'Switch to dark theme';
}
$('theme').onclick = () => {
  const theme = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark';
  document.documentElement.dataset.theme = theme;
  try { localStorage.setItem('scientific-theme', theme); } catch {}
  updateThemeButton();
};
updateThemeButton();
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
    const scroll=element('div','','molecule-scroll');scroll.tabIndex=0;scroll.setAttribute('role','region');scroll.setAttribute('aria-label','Molecular structure; scroll for larger drawings');scroll.append(image);
    const download=element('a','Download SVG');download.href=structure.image_url;download.download='structure-'+(index+1)+'.svg';caption.append(download);figure.append(scroll,caption);gallery.append(figure);
  });container.append(gallery);
}
const basisLabels = {reported:'Reported', computed:'Computed', proposed:'Proposed', unknown:'Unknown', input:'User input'};
function showStructured(view, key = 'science') {
  const section = element('section', '', 'scientific-view');
  if (view.error) { section.append(element('p', view.error, 'scientific-note')); return section; }
  const sources = new Map(view.sources.map(source => [source.id, source]));
  const molecules = new Map(view.molecules.map(molecule => [molecule.id, molecule]));
  const steps = new Map(view.steps.map(step => [step.id, step]));
  function badge(item) { return element('span', basisLabels[item.basis] || item.basis, 'basis basis-' + item.basis); }
  function notes(container, texts) {
    if (!texts.length) return;
    const list = element('ul', '', 'scientific-notes');
    [...new Set(texts)].forEach(text => list.append(element('li', text)));
    container.append(list);
  }
  function links(container, ids) {
    const citations = element('div', '', 'evidence');
    for (const id of new Set(ids)) {
      const source = sources.get(id); if (!source) continue;
      const link = element('a', source.title);
      link.href = source.url || source.artifact_url; link.target = '_blank'; link.rel = 'noopener noreferrer';
      link.title = source.locator; citations.append(link);
    }
    if (citations.children.length) container.append(citations);
  }
  function claim(item) {
    const box = element('div', '', 'attributed-claim');
    box.append(badge(item), element('span', item.text));
    links(box, item.source_ids); notes(box, item.limitations); return box;
  }
  function stepCard(step, routeKey, number) {
    const card = element('article', '', 'scientific-step');
    const heading = element('div', '', 'step-heading');
    heading.append(element('span', String(number), 'step-number'), element('h4', step.title), badge(step));
    card.append(heading);
    if (step.image_url) {
      const figure = element('figure', '', 'scheme-figure');
      const scroll = element('div', '', 'scheme-scroll');
      scroll.tabIndex = 0; scroll.setAttribute('role', 'region'); scroll.setAttribute('aria-label', 'Reaction scheme: ' + step.title);
      const image = element('img', '', 'reaction-scheme'); image.src = step.image_url;
      if (step.scheme_width) image.style.width = step.scheme_width * 0.6 + 'px';
      image.alt = step.reactant_ids.map(id => molecules.get(id).name).join(' + ') + ' → ' + step.product_ids.map(id => molecules.get(id).name).join(' + ') + '. Conditions and yield in step details.';
      image.loading = 'lazy'; scroll.append(image); figure.append(scroll);
      const actions = element('figcaption', '', 'scheme-actions');
      const download = element('a', 'Download SVG'); download.href = step.image_url; download.download = step.id + '-reaction.svg';
      actions.append(download); figure.append(actions); card.append(figure);
    } else card.append(element('p', 'Scheme unavailable for the supplied notation. Structures and conditions are retained in step details.', 'scientific-note'));
    // Keep the scheme compact; the disclosure advertises retained cautions.
    const cautions = new Set([...step.limitations, ...step.conditions.flatMap(item => item.limitations), ...(step.yield_info?.limitations || [])]);
    const detail = element('div', '', 'step-detail');
    notes(detail, step.limitations);
    detail.append(element('p', step.reactant_ids.map(id => molecules.get(id).name).join(' + ') + ' → ' + step.product_ids.map(id => molecules.get(id).name).join(' + ')));
    detail.append(element('h4', 'Conditions'));
    if (step.conditions.length) step.conditions.forEach(item => detail.append(claim(item)));
    else detail.append(element('p', 'Not supplied', 'muted'));
    detail.append(element('h4', 'Yield'));
    detail.append(step.yield_info ? claim(step.yield_info) : element('p', 'Not supplied', 'muted'));
    links(detail, step.source_ids);
    const detailLabel = 'Step details & evidence' + (cautions.size ? ' · ' + cautions.size + (cautions.size === 1 ? ' caution' : ' cautions') : '');
    card.append(disclosure(detailLabel, routeKey + ':' + step.id, detail));
    return card;
  }
  const routed = new Set();
  for (const [index, route] of view.routes.entries()) {
    const box = element('section', '', 'route-section');
    const routeKey = key + ':' + route.id;
    if (route.unreached_target_ids.length) box.append(element('p', 'Incomplete route: does not reach ' + route.unreached_target_ids.map(id => molecules.get(id).name).join(', ') + '.', 'scientific-note'));
    route.step_ids.forEach((id, position) => { box.append(stepCard(steps.get(id), routeKey, position + 1)); routed.add(id); });
    if (view.routes.length === 1) { section.append(element('h3', route.title), box); }
    else {
      const alternative = disclosure(route.title + ' · ' + route.step_ids.length + ' steps', routeKey, box);
      alternative.className += ' route-choice'; alternative.open = index === 0; section.append(alternative);
    }
  }
  const standalone = view.steps.filter(step => !routed.has(step.id));
  if (standalone.length) {
    section.append(element('h3', view.routes.length ? 'Other reaction steps' : 'Reaction schemes'));
    standalone.forEach((step, index) => section.append(stepCard(step, key, index + 1)));
  }
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
    link(row, source.title, source.url || source.artifact_url);
    row.append(element('p', source.locator, 'muted'));
    if (source.url) link(row, 'Saved excerpt', source.artifact_url);
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

function answerDetails(turn) {
  const answer = turn.answer, view = turn.structured_presentation;
  const content = element('div', '', 'answer-detail-content');
  const notes = new Map();
  function addNote(text, scope = '') {
    if (typeof text !== 'string' || !text.trim()) return;
    const key = text.trim().replace(/\s+/g, ' ');
    if (!notes.has(key)) notes.set(key, {text:key, scopes:new Set(), global:false});
    const note = notes.get(key);
    if (scope) note.scopes.add(scope); else note.global = true;
  }
  (answer.uncertainties || []).forEach(text => addNote(text));
  (view?.routes || []).forEach(route => route.limitations.forEach(text => addNote(text, route.title)));
  (view?.molecules || []).forEach(molecule => molecule.limitations.forEach(text => addNote(text, molecule.name)));
  if (notes.size) {
    content.append(element('h4', 'Notes'));
    const list = element('ul', '', 'scientific-notes');
    for (const note of notes.values()) {
      const label = note.global ? '' : [...note.scopes].join('; ') + ': ';
      list.append(element('li', label + note.text));
    }
    content.append(list);
  }
  if (view?.claims?.length) {
    content.append(element('h4', 'Findings'));
    const sources = new Map((view.sources || []).map(source => [source.id, source]));
    for (const claim of view.claims) {
      const row = element('div', '', 'attributed-claim');
      row.append(element('span', basisLabels[claim.basis] || claim.basis, 'basis basis-' + claim.basis), element('span', claim.text));
      for (const text of new Set(claim.limitations)) row.append(element('p', text, 'muted'));
      for (const id of new Set(claim.source_ids)) {
        const source = sources.get(id); if (!source) continue;
        const link = element('a', source.title);
        link.href = source.url || source.artifact_url; link.target = '_blank'; link.rel = 'noopener noreferrer';
        link.title = source.locator; row.append(link);
      }
      content.append(row);
    }
  }
  if (answer.evidence_refs?.length || view?.sources?.length) {
    content.append(element('h4', 'Sources'), sourcesPanel(turn));
  }
  if (!content.children.length) return null;
  const details = disclosure('Notes & sources', turn.id + ':details', content);
  details.className += ' answer-details';
  return details;
}

function answerCard(turn) {
  const card = element('article', '', 'assistant');
  if (turn.progress?.length || runningStates.has(turn.status)) card.append(investigationTimeline(turn));
  if (turn.answer) {
    const answer = turn.answer;
    const view = turn.structured_presentation;
    const routeFirst = !view?.error && view?.steps?.length > 0;
    if (!routeFirst) card.append(formattedMessage(answer.answer_markdown, turn.answer_presentation));
    if (view && (view.error || view.steps?.length)) {
      card.append(showStructured(view, turn.id + ':science'));
    }
    if (routeFirst) card.append(formattedMessage(answer.answer_markdown, turn.answer_presentation));
    const details = answerDetails(turn);
    if (details) card.append(details);
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
    const caption = element('p', '', 'progress-caption'); caption.hidden = true;
    caption.id = 'progress-caption';
    progress.append(line, caption);
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
    const position = $('scroll-area').scrollTop;
    const expanded = new Map(Array.from($('messages').querySelectorAll('details'))
      .filter(node => node.dataset.key).map(node => [node.dataset.key, Boolean(node.open)]));
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
    $('messages').querySelectorAll('details').forEach(node => {
      if (expanded.has(node.dataset.key)) node.open = expanded.get(node.dataset.key);
    });
    displayed = signature;
    if (follow) requestAnimationFrame(toBottom);
    else $('scroll-area').scrollTop = position;
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
  const names = {command_execution: 'Command', web_search: 'Web search', mcp_tool_call: 'Scientific tool',
    file_change: 'Saving investigation files', todo_list: 'Updating the plan'};
  const label = event.title || names[event.kind] || 'Agent activity';
  const failed = ['failed', 'item.failed', 'error', 'timed_out'].includes(event.status) || (Number.isInteger(event.exit_code) && event.exit_code !== 0);
  const status = failed ? 'failed' : ['completed', 'item.completed'].includes(event.status) ? 'finished' :
    ['cancelled', 'canceled'].includes(event.status) ? 'stopped' : 'running';
  return label + ' · ' + status + (Number.isInteger(event.exit_code) ? ' · exit ' + event.exit_code : '');
}
function investigationTimeline(turn) {
  const events = element('ol', '', 'investigation-timeline');
  events.dataset.turn = turn.id;
  events.setAttribute('aria-label', 'Investigation progress');
  if (runningStates.has(turn.status)) {
    events.id = 'live-timeline';
    events.setAttribute('role', 'log'); events.setAttribute('aria-live', 'polite');
    events.setAttribute('aria-relevant', 'additions text');
  }
  renderTimeline(events, turn.progress || []);
  return events;
}
function actionSummary(event) {
  const subject = (event.detail || event.title || '').replace(/\s+/g, ' ').trim();
  const failed = ['failed', 'item.failed', 'error', 'timed_out'].includes(event.status) ||
    (Number.isInteger(event.exit_code) && event.exit_code !== 0);
  if (failed) return 'Failed: ' + (subject || 'tool action');
  if (['cancelled', 'canceled'].includes(event.status)) return 'Stopped: ' + (subject || 'tool action');
  const finished = ['completed', 'item.completed'].includes(event.status);
  const verbs = {
    command_execution: ['Running', 'Ran'], custom_execution: ['Running', 'Ran'],
    web_search: ['Searching', 'Searched'], file_change: ['Updating', 'Updated'],
    mcp_tool_call: ['Calling', 'Called'], scientific_call: ['Checking', 'Checked'],
    literature_source: ['Reading', 'Read'], todo_list: ['Planning', 'Planned'],
  };
  return (verbs[event.kind] || ['Working on', 'Finished'])[finished ? 1 : 0] + ' ' + (subject || 'investigation');
}
function setTimelineText(node, text) {
  if (node.textContent !== text) node.textContent = text;
}
function createTimelineRow(key, update) {
  const row = element('li', '', 'timeline-entry');
  row.dataset.key = key;
  row.dataset.mode = update ? 'update' : 'action';
  if (update) {
    row.append(element('p', '', 'investigation-update'));
  } else {
    const action = element('details', '', 'timeline-action');
    action.dataset.key = key; action.open = false;
    const summary = element('summary');
    const icon = element('span', '', 'tool-icon'); icon.setAttribute('aria-hidden', 'true');
    summary.append(icon, element('span', '', 'tool-summary'));
    const content = element('div', '', 'action-details');
    const metadata = element('div', '', 'action-metadata');
    metadata.append(element('span', '', 'activity-status'));
    content.append(metadata, element('div', '', 'activity-detail'), element('div', '', 'activity-error'));
    action.append(summary, content);
    row.append(action);
  }
  return row;
}
function updateTimelineRow(row, event) {
  const signature = JSON.stringify(event);
  if (row.dataset.signature === signature) return;
  if (event.kind === 'agent_update') {
    setTimelineText(row.firstElementChild, event.detail || '');
  } else {
    const action = row.firstElementChild;
    const summary = action.firstElementChild;
    const content = action.lastElementChild;
    const metadata = content.firstElementChild;
    const icons = {command_execution:'›_', custom_execution:'›_', web_search:'⌕', file_change:'✎',
      mcp_tool_call:'◇', scientific_call:'◇', literature_source:'↗', todo_list:'☷'};
    const label = actionSummary(event);
    action.dataset.failed = String(label.startsWith('Failed:'));
    summary.title = label;
    setTimelineText(summary.firstElementChild, icons[event.kind] || '·');
    setTimelineText(summary.lastElementChild, label);
    setTimelineText(metadata.firstElementChild, activityLabel(event));
    const at = event.at ? new Date(event.at) : null;
    let timestamp = metadata.children[1];
    if (at && !Number.isNaN(at.getTime())) {
      if (!timestamp) { timestamp = element('time'); metadata.append(timestamp); }
      setTimelineText(timestamp, at.toLocaleTimeString([], {hour:'2-digit', minute:'2-digit', second:'2-digit'}));
      timestamp.setAttribute('datetime', at.toISOString());
    } else if (timestamp) timestamp.remove();
    const detail = event.detail || (!event.title ? 'The runtime did not provide action details.' : '');
    setTimelineText(content.children[1], detail);
    content.children[1].hidden = !detail;
    setTimelineText(content.children[2], event.failure_detail || '');
    content.children[2].hidden = !event.failure_detail;
  }
  row.dataset.signature = signature;
}
function renderTimeline(events, progress) {
  const signature = JSON.stringify(progress);
  if (events.dataset.signature === signature) return false;
  const rows = new Map(Array.from(events.children).map(row => [row.dataset.key, row]));
  const retained = new Set();
  for (const [index, event] of progress.entries()) {
    const key = events.dataset.turn + ':action:' + (event.activity_id || event.item_id || index);
    const update = event.kind === 'agent_update';
    let row = rows.get(key);
    if (row && row.dataset.mode !== (update ? 'update' : 'action')) { row.remove(); row = null; }
    if (!row) row = createTimelineRow(key, update);
    updateTimelineRow(row, event);
    retained.add(row);
    // Preserve disclosure focus, selection, and live-region history during polling.
    if (events.children[index] !== row) events.insertBefore(row, events.children[index] || null);
  }
  for (const row of Array.from(events.children)) if (!retained.has(row)) row.remove();
  events.dataset.signature = signature;
  return true;
}
function updateProgress() {
  const turn = active?.conversation_id === identity ? active : conversation?.turns.at(-1);
  const label = $('progress-label');
  if (label && turn && runningStates.has(turn.status)) {
    label.textContent = stopping === identity ? 'Stopping investigation…' :
      turn.status === 'preparing' ? 'Preparing your workspace…' : turn.status === 'queued' ? 'Starting investigation…' :
      turn.repair_attempts ? 'Refining the answer…' : 'Investigating…';
    $('elapsed').textContent = elapsedText(turn.created_at);
    const caption = $('progress-caption');
    if (caption) {
      const times = (turn.progress || []).map(event => event.updated_at || event.at).filter(value => Number.isFinite(Date.parse(value)));
      const latest = times.sort((a, b) => Date.parse(a) - Date.parse(b)).at(-1);
      const since = latest || turn.created_at;
      const quiet = turn.status === 'running' && Date.now() - Date.parse(since) >= 60000;
      caption.hidden = !quiet;
      caption.textContent = quiet ?
        (latest ? 'Last recorded activity ' : 'No recorded activity for ') + elapsedText(since) +
          (latest ? ' ago. ' : '. ') + 'The runtime is still running; no new update has arrived.' :
        '';
    }
    const events = $('live-timeline');
    const follow = nearBottom();
    const position = $('scroll-area').scrollTop;
    if (events && renderTimeline(events, turn.progress || [])) {
      if (follow) requestAnimationFrame(toBottom);
      else $('scroll-area').scrollTop = position;
    }
  }
  if (active) $('active-text').textContent = (stopping ? 'Stopping investigation' : 'An investigation is running') + ' · ' + elapsedText(active.created_at);
}
function resizeInput() {
  $('question').style.height = 'auto'; $('question').style.height = Math.min($('question').scrollHeight, 140) + 'px';
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
    if (config.research_settings) {
      const settings = config.research_settings;
      $('environment').append(element('p', 'Research settings: ' + settings.profile +
        ' · reasoning ' + (settings.requested.reasoning_effort || 'inherited') +
        ' · web search ' + (settings.requested.web_search || 'inherited') +
        ' (requested; access is checked during each investigation).'));
    }
    active = (await request('/activity')).active; await list(); loaded = true; updateControls();
    const id = location.hash.slice(1);
    if (/^[0-9a-f]{32}$/.test(id)) await openConversation(id);
    else if (active) await openConversation(active.conversation_id);
    else await poll();
  } catch (error) { $('error').textContent = 'Unable to connect. ' + error.message; updateControls(); }
}
setInterval(updateProgress, 1000);
initialize();
