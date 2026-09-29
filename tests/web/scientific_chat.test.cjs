// Controller tests with a small DOM double: no browser, network, or model calls.
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '../..');

class Node {
  constructor(tag = 'div') {
    this.tag = tag; this.children = []; this.parentNode = null; this.dataset = {}; this.style = {};
    this.attributes = {}; this.value = ''; this.hidden = false; this.disabled = false;
    this.classList = {toggle() {}}; this.scrollHeight = 400; this.scrollTop = 0; this.clientHeight = 400;
  }
  set textContent(value) { this.replaceChildren(); this.text = value; }
  get textContent() { return this.text || ''; }
  get firstElementChild() { return this.children[0] || null; }
  get lastElementChild() { return this.children.at(-1) || null; }
  append(...nodes) { for (const node of nodes) this.insertBefore(node, null); }
  replaceChildren(...nodes) {
    for (const child of [...this.children]) child.remove();
    this.append(...nodes);
  }
  insertBefore(node, reference) {
    if (node === reference) return node;
    if (reference !== null && !this.children.includes(reference)) throw new Error('Reference is not a child');
    node.remove();
    this.children.splice(reference === null ? this.children.length : this.children.indexOf(reference), 0, node);
    node.parentNode = this;
    return node;
  }
  remove() {
    if (this.parentNode) {
      const siblings = this.parentNode.children;
      siblings.splice(siblings.indexOf(this), 1);
      this.parentNode = null;
    }
  }
  setAttribute(name, value) { this.attributes[name] = value; }
  focus() {}
  requestSubmit() { this.submissions = (this.submissions || 0) + 1; }
  descendants() { return this.children.flatMap(child => [child, ...child.descendants()]); }
  querySelectorAll(selector) {
    return this.descendants().filter(node => selector === 'details[open]' ? node.tag === 'details' && node.open : node.tag === selector);
  }
}

function harness() {
  const html = fs.readFileSync(path.join(root, 'app/web_api/scientific_chat.html'), 'utf8');
  const ids = Object.fromEntries([...html.matchAll(/id="([^"]+)"/g)].map(match => [match[1], new Node()]));
  const context = vm.createContext({
    document: {
      documentElement: {dataset: {theme: 'dark'}},
      getElementById: id => ids[id] || Object.values(ids).flatMap(node => node.descendants()).find(node => node.id === id) || null,
      createElement: tag => new Node(tag), createTextNode: text => { const node = new Node(); node.textContent = text; return node; },
      querySelectorAll: () => [], addEventListener() {},
    },
    location: {hash: ''}, window: {addEventListener() {}}, navigator: {clipboard: {writeText: async () => {}}},
    matchMedia: () => ({matches: false}), requestAnimationFrame: callback => callback(),
    setTimeout: () => 1, clearTimeout() {}, setInterval() {},
  });
  const script = fs.readFileSync(path.join(root, 'app/web_api/scientific_chat.js'), 'utf8');
  vm.runInContext(script.replace(/initialize\(\);\s*$/, ''), context);
  function run(code) { return vm.runInContext(code, context); }
  function respond(handler) { context.respond = handler; run('request = (path, options) => respond(path, options)'); }
  return {ids, context, run, respond};
}

test('Stop remains available in another chat and cancels the actual owner', async () => {
  const {ids, run, respond} = harness();
  run("loaded=true; identity='history'; active={conversation_id:'working',status:'running'}; conversation={turns:[]};");
  ids.question.value = 'A follow-up draft';
  run('updateControls()');
  assert.equal(ids.send.hidden, true);
  assert.equal(ids.cancel.hidden, false);
  assert.equal(ids['global-activity'].hidden, false);
  assert.equal(ids.question.disabled, false);
  const calls = [];
  respond(async (route, options) => {
    calls.push([route, options?.method]);
    if (route.endsWith('/cancel')) return {cancellation_requested: true};
    if (route === '/activity') return {active: null};
    if (route === '/conversations') return [];
    throw new Error('Unexpected route: ' + route);
  });
  await run('stopActive()');
  assert.deepEqual(calls[0], ['/conversations/working/cancel', 'POST']);
  assert.equal(ids.cancel.hidden, true);
  assert.equal(ids.send.disabled, false);
  assert.equal(ids.question.value, 'A follow-up draft');
});

function byClass(node, name) {
  return node.descendants().filter(child => (child.className || '').split(/\s+/).includes(name));
}

test('one timeline retains older entries and the main conversation preserves reading position', () => {
  const {run, context, ids} = harness();
  context.progress = Array.from({length: 45}, (_, index) => ({
    activity_id:'action:' + index, kind:'command_execution', status:'completed',
    at:'2026-09-27T13:16:59Z', detail:'python analysis_' + index + '.py', exit_code:0,
  }));
  run("identity='chat'; active={conversation_id:'chat',status:'running',progress}; renderConversation({title:'Work',turns:[{id:'t',status:'running',question:'Analyze',progress}]})");
  const timeline = run("$('live-timeline')");
  assert.equal(timeline.children.length, 45);
  assert.ok(timeline.children[0].descendants().some(node => node.textContent === 'python analysis_0.py'));
  ids['scroll-area'].scrollHeight = 2200; ids['scroll-area'].clientHeight = 220;
  ids['scroll-area'].scrollTop = 80;
  context.progress.push({activity_id:'search',kind:'web_search',status:'in_progress',detail:'Suzuki reaction conditions'});
  run('updateProgress()');
  assert.equal(timeline.children.length, 46);
  assert.equal(ids['scroll-area'].scrollTop, 80);
  const firstRow = timeline.children[0];
  run('updateProgress()');
  assert.equal(timeline.children[0], firstRow, 'unchanged polling must not rebuild the timeline');
  ids['scroll-area'].scrollTop = 1980;
  context.progress.push({activity_id:'search-2',kind:'web_search',status:'completed',detail:'Published preparation'});
  run('updateProgress()');
  assert.equal(ids['scroll-area'].scrollTop, ids['scroll-area'].scrollHeight);
});

test('completed timeline keeps command failures literal without an activity box or debug download', () => {
  const {run} = harness();
  const card = run(`answerCard({id:'done',status:'completed',debug_log_available:true,answer:{answer_markdown:'Result'},
    progress:[{activity_id:'cmd',kind:'command_execution',status:'completed',exit_code:1,detail:'<script>alert(1)</script>'}]})`);
  const timeline = card.children[0];
  assert.equal(timeline.className, 'investigation-timeline');
  const action = byClass(timeline, 'timeline-action')[0];
  assert.equal(action.dataset.key, 'done:action:cmd');
  assert.equal(Boolean(action.open), false);
  assert.match(byClass(action, 'tool-summary')[0].textContent, /failed/i);
  assert.equal(byClass(action, 'activity-detail')[0].textContent, '<script>alert(1)</script>');
  assert.equal(card.querySelectorAll('script').length, 0);
  assert.ok(card.children.some(node => node.textContent === 'Result'));
  assert.ok(!card.querySelectorAll('summary').some(node => /View activity|Investigation updates/.test(node.textContent)));
  assert.ok(!card.querySelectorAll('a').some(node => node.textContent === 'Download debug log'));
});

test('inline action details keep diagnostics and omit unknown timestamps', () => {
  const {run, context} = harness();
  context.timeline = new Node('ol'); context.timeline.dataset.turn = 't';
  run(`renderTimeline(timeline, [{activity_id:'read',kind:'command_execution',title:'Read old.md',status:'failed',exit_code:1,
    at:null,detail:'Get-Content old.md',failure_detail:'Cannot find path <old.md>'}])`);
  const action = byClass(context.timeline, 'timeline-action')[0];
  assert.equal(action.querySelectorAll('time').length, 0);
  assert.equal(byClass(action, 'activity-status')[0].textContent, 'Read old.md · failed · exit 1');
  assert.equal(byClass(action, 'activity-error')[0].textContent, 'Cannot find path <old.md>');
  assert.equal(action.querySelectorAll('old.md').length, 0);
});

test('running commands use one readable summary line while retaining expanded command text', () => {
  const {run, context} = harness();
  context.command = 'python -c "\n  import json\n  print(json.dumps({\"ok\": True}))\n"';
  context.timeline = new Node('ol'); context.timeline.dataset.turn = 't';
  run("renderTimeline(timeline, [{activity_id:'cmd',kind:'command_execution',status:'in_progress',detail:command}])");
  const action = byClass(context.timeline, 'timeline-action')[0];
  const summary = byClass(action, 'tool-summary')[0].textContent;
  assert.match(summary, /^Running python -c/);
  assert.doesNotMatch(summary, /\n|\r|\s{2}/);
  assert.equal(byClass(action, 'activity-detail')[0].textContent, context.command);
  assert.equal(action.open, false);
});

test('public prose and collapsed actions are interleaved in their recorded order', () => {
  const {run, context} = harness();
  context.turn = {id:'t',status:'running',debug_log_available:true,progress:[
    {activity_id:'note-1',kind:'agent_update',detail:'I found a preparation; I’m checking its yield.'},
    {activity_id:'search',kind:'web_search',title:'Search the web: patent Example 2',status:'in_progress',detail:'patent Example 2'},
    {activity_id:'note-2',kind:'agent_update',detail:'The reported yield applies to a different substrate.'},
  ]};
  const card = run('answerCard(turn)');
  const timeline = byClass(card, 'investigation-timeline')[0];
  assert.equal(timeline.children.length, 3);
  assert.equal(byClass(timeline.children[0], 'investigation-update')[0].textContent, context.turn.progress[0].detail);
  assert.equal(byClass(timeline.children[1], 'timeline-action').length, 1);
  assert.equal(byClass(timeline.children[2], 'investigation-update')[0].textContent, context.turn.progress[2].detail);
  context.turn.status = 'failed';
  const saved = run('answerCard(turn)');
  const savedTimeline = byClass(saved, 'investigation-timeline')[0];
  assert.equal(savedTimeline.children.length, 3);
  assert.equal(byClass(savedTimeline, 'investigation-update').length, 2);
  assert.ok(!saved.querySelectorAll('details').some(node => /:(activity|updates)$/.test(node.dataset.key || '')));
});

test('lifecycle and new events retain focused disclosure nodes and expanded diagnostics', () => {
  const {run, context} = harness();
  context.timeline = new Node('ol'); context.timeline.dataset.turn = 't';
  context.progress = [
    {activity_id:'intro',kind:'agent_update',detail:'I’ll inspect the target.'},
    {activity_id:'cmd',kind:'command_execution',title:'Run Python script: inspect.py',status:'in_progress',detail:'python inspect.py'},
    {activity_id:'next',kind:'agent_update',detail:'The target identity is clear.'},
  ];
  run('renderTimeline(timeline, progress)');
  const originalRow = context.timeline.children[1];
  const originalAction = byClass(originalRow, 'timeline-action')[0];
  const originalSummary = originalAction.firstElementChild;
  originalAction.open = true;
  Object.assign(context.progress[1], {status:'completed',exit_code:0});
  run('renderTimeline(timeline, progress)');
  assert.equal(context.timeline.children.length, 3);
  const action = byClass(context.timeline.children[1], 'timeline-action')[0];
  assert.equal(context.timeline.children[1], originalRow);
  assert.equal(action, originalAction);
  assert.equal(action.firstElementChild, originalSummary, 'status updates must retain the keyboard-focus target');
  assert.equal(action.open, true);
  assert.match(byClass(action, 'activity-status')[0].textContent, /finished/);
  assert.equal(byClass(context.timeline.children[0], 'investigation-update')[0].textContent, context.progress[0].detail);
  const first = context.timeline.children[0];
  run('renderTimeline(timeline, progress)');
  assert.equal(context.timeline.children[0], first);
  context.progress.push({activity_id:'new',kind:'web_search',status:'in_progress',detail:'Read a primary source'});
  run('renderTimeline(timeline, progress)');
  assert.equal(byClass(context.timeline, 'timeline-action')[0].firstElementChild, originalSummary);
  assert.equal(originalAction.open, true);
  // Reconciliation moves an existing row without replacing its controls.
  context.progress = [context.progress[1], context.progress[2], context.progress[3]];
  run('renderTimeline(timeline, progress)');
  assert.equal(context.timeline.children.length, 3);
  assert.equal(context.timeline.children[0], originalRow);
  assert.equal(originalAction.firstElementChild, originalSummary);
  assert.equal(originalAction.open, true);
  assert.equal(first.parentNode, null, 'removed history entries leave the DOM');
});

test('finishing a turn preserves expanded actions and the main conversation reading position', () => {
  const {run, ids} = harness();
  run(`identity='chat'; renderConversation({title:'Activity',turns:[{id:'turn',status:'running',
    question:'Analyze',progress:[{activity_id:'search',kind:'web_search',status:'completed',detail:'Reaction conditions'}]}]})`);
  const action = byClass(ids.messages, 'timeline-action')[0];
  action.open = true;
  ids['scroll-area'].scrollHeight = 2200; ids['scroll-area'].clientHeight = 220;
  ids['scroll-area'].scrollTop = 90;
  run(`renderConversation({title:'Activity',turns:[{id:'turn',status:'completed',question:'Analyze',
    answer:{answer_markdown:'Done'},progress:[{activity_id:'search',kind:'web_search',status:'completed',detail:'Reaction conditions'}]}]})`);
  assert.equal(byClass(ids.messages, 'timeline-action')[0].open, true);
  assert.equal(ids['scroll-area'].scrollTop, 90);
  assert.equal(byClass(ids.messages, 'investigation-timeline').length, 1);
});

test('a rejected cancellation restores Stop and displays the error', async () => {
  const {ids, run, respond} = harness();
  run("loaded=true; active={conversation_id:'working'};");
  respond(async () => { throw new Error('Reload to obtain the session token'); });
  await run('stopActive()');
  assert.equal(ids.cancel.disabled, false);
  assert.equal(ids.cancel.hidden, false);
  assert.match(ids.error.textContent, /Reload/);
});

test('late navigation responses cannot replace the selected conversation', async () => {
  const {run, respond} = harness();
  let resolve;
  respond(() => new Promise(done => { resolve = done; }));
  run("identity='first'; navigation=1");
  const loading = run("loadConversation('first', 1)");
  run("identity='second'; navigation=2; conversation={title:'Second',turns:[]}");
  resolve({title:'First', turns:[]});
  await loading;
  assert.equal(run('conversation.title'), 'Second');
});

test('Enter sends, Shift+Enter and IME composition keep editing', () => {
  const {ids} = harness();
  ids.send.disabled = false;
  let prevented = 0;
  const event = {key:'Enter', preventDefault:() => prevented++};
  ids.question.onkeydown({...event, shiftKey:true});
  ids.question.onkeydown({...event, isComposing:true});
  assert.equal(prevented, 0);
  ids.question.onkeydown(event);
  assert.equal(ids.form.submissions, 1);
  ids.send.disabled = true;
  ids.question.onkeydown(event);
  assert.equal(ids.form.submissions, 1);
});

test('new chat preserves drafts and global running activity', () => {
  const {ids, run, respond} = harness();
  respond(async () => []);
  run("loaded=true; identity='old'; active={conversation_id:'working'};");
  ids.question.value = 'Unsaved draft';
  run('newChat()');
  assert.equal(run("drafts.get('old')"), 'Unsaved draft');
  assert.equal(ids.question.value, '');
  assert.equal(ids['global-activity'].hidden, false);
  assert.equal(ids.cancel.hidden, false);
});

test('legacy answers have one compact notes and sources footer without molecule galleries', () => {
  const {run} = harness();
  const card = run(`answerCard({id:'one',status:'completed',answer:{answer_markdown:'A result',evidence_refs:['sha256:abc'],uncertainties:['Unverified']},
    answer_presentation:{structures:[{smiles:'CCO',image_url:'data:image/svg+xml;base64,abc'}]}})`);
  assert.equal(card.children[0].textContent, 'A result');
  assert.equal(card.querySelectorAll('details').length, 1);
  assert.equal(card.querySelectorAll('details[open]').length, 0);
  assert.equal(card.querySelectorAll('summary')[0].textContent, 'Notes & sources');
  assert.ok(card.descendants().some(node => node.textContent === 'Unverified'));
  assert.ok(card.querySelectorAll('a').some(node => node.textContent === 'Recorded evidence 1'));
  assert.equal(card.querySelectorAll('img').length, 0);
});

test('supporting reactions retain source identity, separate observations and evidence gaps', () => {
  const {run, context} = harness();
  const step = {id:'s1',title:'Proposed step',basis:'proposed',reactant_ids:['a'],product_ids:['b'],
    after_step_ids:[],conditions:[],yield_info:null,source_ids:[],limitations:[]};
  context.supportView = {sources:[],molecules:[{id:'a',name:'A'},{id:'b',name:'B'}],routes:[],steps:[step]};
  let card = run('showStructured(supportView)');
  assert.ok(card.descendants().some(node => node.textContent === 'No supporting-reaction inspection attached to this step.'));
  step.supporting_evidence = [{artifact_ref:'sha256:fixture',artifact_url:'/saved/inspection',status:'precedents_available',
    scope:'selected_template_top_20',saved_match_count:4,distinct_references_on_page:1,page:{next_offset:1},
    precedents:[{match_id:'match1',reaction_id:'reaction-1',reference_title:'Reported patent example',reference_url:'https://example.org/patent',
      reaction_smiles:'CCO>>CC=O',image_url:'data:image/svg+xml;base64,abc',same_recorded_product:false,same_recorded_precursors:false,
      product_similarity:0.4,precursor_similarity:0.3,limitations:['Transfer remains uncertain.'],
      observations:[{observation_id:'obs1',yield_pct:40,resolved_recipe:{solvents:[{canonical_name:'Ethanol'}]}},{observation_id:'obs2',yield_pct:null}],
      procedures:[{observation_id:'obs1',procedure_text:'Literal source <script>text</script>'}]}]}];
  card = run('showStructured(supportView)');
  assert.equal(card.querySelectorAll('img').length,1);
  assert.ok(card.querySelectorAll('summary').some(node => node.textContent === 'Supporting reactions (1)'));
  assert.ok(card.querySelectorAll('a').some(node => node.href === 'https://example.org/patent'));
  assert.ok(card.descendants().some(node => node.textContent === 'Reported yield: 40%'));
  assert.ok(card.descendants().some(node => node.textContent === 'Reported yield: not supplied'));
  assert.ok(card.querySelectorAll('code').some(node => node.textContent === 'CCO>>CC=O'));
  const conditions = card.descendants().find(node => node.className === 'precedent-conditions');
  assert.equal(conditions.parentNode.className, 'precedent-card');
  assert.ok(conditions.descendants().some(node => node.textContent === 'Reagents & solvents: Ethanol'));
  const cautions = card.querySelectorAll('summary').find(node => node.textContent === 'Match details & cautions');
  assert.ok(!cautions.parentNode.open);
  assert.ok(card.descendants().some(node => node.textContent === 'Literal source <script>text</script>'));
  assert.ok(!card.querySelectorAll('a').some(node => node.download));
  step.supporting_evidence = [{status:'no_precedents_retrieved',precedents:[],scope:'saved_assessment_matches',saved_match_count:0,distinct_references_on_page:0}];
  card = run('showStructured(supportView)');
  assert.ok(card.descendants().some(node => node.textContent === 'No supporting experimental precedent retrieved in this inspection.'));
});

test('reaction schemes are visible with details collapsed and references human-readable', () => {
  const {run, context} = harness();
  const attribution = {basis:'proposed', source_ids:['paper'], limitations:['Feasibility unverified']};
  context.fixture = {
    id:'scheme', answer:{answer_markdown:'A concise conclusion', evidence_refs:[], uncertainties:[]},
    structured_presentation:{
      sources:[{id:'paper',title:'Patent Example 2',url:'https://example.org/patent',artifact_url:'/saved/excerpt',locator:'Example 2'}],
      molecules:[{id:'a',name:'Reactant',smiles:'CCO',...attribution},{id:'b',name:'Product',smiles:'CC=O',...attribution}],
      target_molecule_ids:['b'], routes:[], claims:[],
      steps:[{id:'s1',title:'Oxidation hypothesis',reactant_ids:['a'],product_ids:['b'],after_step_ids:[],
        conditions:[{text:'Unknown oxidant',...attribution}],yield_info:null,reaction_smiles:'CCO>>CC=O',
        image_url:'data:image/svg+xml;base64,abc',scheme_width:980,...attribution}],
    },
  };
  const card = run('answerCard(fixture)');
  const view = card.children.find(node => node.className === 'scientific-view');
  assert.ok(view, 'schemes should not be inside a collapsed details wrapper');
  assert.equal(card.children[0], view, 'show the explicit route before its concise explanation');
  assert.equal(card.children[1].textContent, 'A concise conclusion');
  assert.equal(card.querySelectorAll('details[open]').length, 0);
  const images = card.querySelectorAll('img');
  assert.equal(images.length, 1);
  assert.equal(images[0].style.width, '588px', 'prefer a compact preview at 60% of native SVG width');
  assert.match(images[0].alt, /Reactant → Product/);
  assert.equal(images[0].src, context.fixture.structured_presentation.steps[0].image_url);
  assert.equal(card.querySelectorAll('a').some(node => node.download), false,
    'reaction schemes should not show a download button');
  const stepDetails = card.querySelectorAll('details').find(node => node.dataset.key === 'scheme:science:s1');
  assert.equal(stepDetails.firstElementChild.textContent, 'Step details & evidence · 1 caution');
  assert.ok(stepDetails.descendants().some(node => node.textContent === 'Feasibility unverified'));
  assert.equal(byClass(view, 'scientific-step')[0].children.filter(node => node.tag === 'ul').length, 0,
    'detailed cautions should not duplicate the concise explanation below the scheme');
  assert.ok(stepDetails.descendants().some(node => node.textContent === 'Unknown oxidant'));
  assert.ok(stepDetails.descendants().some(node => node.textContent === 'Reactant → Product'));
  assert.ok(byClass(card, 'step-heading')[0].descendants().some(node => node.textContent === 'Proposed'),
    'the actual step attribution remains visible');
  assert.ok(!card.descendants().some(node => node.textContent === 'Synthetic direction' ||
    node.textContent.startsWith('Reaction schemes · declared structures')));
  assert.ok(card.querySelectorAll('a').some(node => node.href === 'https://example.org/patent' && node.textContent === 'Patent Example 2'));
  assert.ok(card.querySelectorAll('a').some(node => node.href === '/saved/excerpt' && node.textContent === 'Saved excerpt'));
  assert.equal(card.querySelectorAll('code').length, 0, 'raw SMILES are retained in data, not duplicated in the answer UI');
  assert.ok(card.descendants().some(node => node.tag === 'li' && node.textContent === 'Feasibility unverified'));
  context.fixture.structured_presentation.molecules.push({id:'c',name:'Final product',smiles:'CC(=O)O',...attribution});
  context.fixture.structured_presentation.steps.push({
    ...context.fixture.structured_presentation.steps[0], id:'s2', title:'Further oxidation',
    reactant_ids:['b'], product_ids:['c'], after_step_ids:['s1'], reaction_smiles:'CC=O>>CC(=O)O',
  });
  context.fixture.structured_presentation.target_molecule_ids = ['c'];
  context.fixture.structured_presentation.routes = [
    {id:'r1',title:'First hypothesis',step_ids:['s1','s2'],limitations:['Supply unconfirmed'],unreached_target_ids:[]},
    {id:'r2',title:'Alternative hypothesis',step_ids:['s1'],limitations:[],unreached_target_ids:['c']},
  ];
  const routedCard = run('answerCard(fixture)');
  const alternatives = routedCard.querySelectorAll('details').filter(node => node.className.includes('route-choice'));
  assert.equal(alternatives.length, 2);
  assert.equal(alternatives[0].open, true);
  assert.equal(alternatives[1].open, false);
  assert.equal(byClass(alternatives[0], 'scientific-step').length, 2, 'keep both steps of a connected route visible');
  const footer = byClass(routedCard, 'answer-details')[0];
  assert.ok(footer.descendants().some(node => node.textContent === 'First hypothesis: Supply unconfirmed'));
  assert.equal(Boolean(footer.open), false);
  assert.equal(routedCard.children.filter(node => node.tag === 'details').length, 1);
  assert.ok(!routedCard.querySelectorAll('summary').some(node => /Step connections|Molecules & SMILES|Route limitations|Uncertainty &|Additional scientific/.test(node.textContent)));
  assert.ok(alternatives[1].descendants().some(node => node.textContent.includes('Incomplete route: does not reach Final product')));
  context.fixture.question = 'Investigate this route';
  context.fixture.status = 'completed';
  run("renderConversation({title:'Routes',turns:[fixture]})");
  const first = run("$('messages').querySelectorAll('details').find(node => node.className.includes('route-choice'))");
  first.open = false;
  context.fixture.answer_ref = 'new-presentation';
  run("renderConversation({title:'Routes',turns:[fixture]})");
  const retained = run("$('messages').querySelectorAll('details').find(node => node.className.includes('route-choice'))");
  assert.equal(retained.open, false, 'a user-collapsed route stays closed when the conversation updates');
});

test('one footer deduplicates notes while preserving scoped caveats, claims and source links', () => {
  const {run, context} = harness();
  context.fixture = {
    id:'notes', answer:{answer_markdown:'Brief explanation', evidence_refs:[], uncertainties:['Stock unknown', '<script>plain text</script>']},
    structured_presentation:{
      sources:[{id:'paper',title:'Example 4',url:'https://example.org/patent',artifact_url:'/saved/source',locator:'Example 4'}],
      molecules:[{id:'a',name:'Precursor',smiles:'CCO',basis:'proposed',limitations:['Stock  unknown'],source_ids:[]}],
      routes:[{id:'r1',title:'Route A',step_ids:[],limitations:['Stock unknown', 'Specific precursor needed'],unreached_target_ids:[]}],
      steps:[], claims:[{text:'Conflicting yield entries',basis:'reported',source_ids:['paper'],limitations:['Yield omitted']}],
    },
  };
  const before = JSON.stringify(context.fixture);
  const card = run('answerCard(fixture)');
  assert.equal(card.children[0].textContent, 'Brief explanation');
  assert.equal(byClass(card, 'scientific-view').length, 0, 'no empty structure section for non-route answers');
  const footer = card.querySelectorAll('details');
  assert.equal(footer.length, 1);
  const list = footer[0].querySelectorAll('li');
  assert.equal(list.filter(node => node.textContent === 'Stock unknown').length, 1);
  assert.ok(list.some(node => node.textContent === 'Route A: Specific precursor needed'));
  assert.ok(list.some(node => node.textContent === '<script>plain text</script>'));
  assert.equal(card.querySelectorAll('script').length, 0);
  assert.ok(footer[0].descendants().some(node => node.textContent === 'Reported'));
  assert.ok(footer[0].descendants().some(node => node.textContent === 'Conflicting yield entries'));
  assert.ok(footer[0].querySelectorAll('a').some(node => node.href === '/saved/source'));
  assert.equal(JSON.stringify(context.fixture), before, 'display grouping must not rewrite saved evidence');
  assert.equal(run("answerCard({id:'empty',answer:{answer_markdown:'Short answer'}})").querySelectorAll('details').length, 0);
});

test('progress uses recorded events, elapsed time, and keeps action details open across updates', () => {
  const {ids, run} = harness();
  run("identity='working'; active={conversation_id:'working', status:'running', created_at:new Date(Date.now()-65000).toISOString(),progress:[{kind:'web_search',status:'item.completed',at:new Date().toISOString()}]}; renderConversation({title:'Question',turns:[{id:'t1',question:'Question',status:'running'}]})");
  const detail = ids.messages.querySelectorAll('details')[0]; detail.open = true;
  assert.match(run("$('elapsed').textContent"), /1m 5s/);
  assert.equal(run("$('progress-label').textContent"), 'Investigating…');
  assert.ok(byClass(ids.messages, 'activity-status').some(node => node.textContent === 'Web search · finished'));
  run("active.progress.push({kind:'command_execution',status:'in_progress'}); renderConversation({title:'Question',turns:[{id:'t1',question:'Question',status:'running'}]})");
  const refreshedAction = ids.messages.querySelectorAll('details')[0];
  assert.equal(refreshedAction.dataset.key, detail.dataset.key);
  assert.equal(refreshedAction.open, true);
});

test('quiet runs report the last actual activity without manufacturing agent commentary', () => {
  const {run} = harness();
  run("identity='working'; active={conversation_id:'working',status:'running',created_at:new Date(Date.now()-180000).toISOString(),progress:[{kind:'scientific_call',status:'failed',title:'Scientific call: assess_route_proposal',at:new Date(Date.now()-65000).toISOString(),failure_detail:'Invalid evidence reference'}]}; renderConversation({title:'Question',turns:[{id:'t',question:'Question',status:'running'}]})");
  assert.match(run("$('progress-caption').textContent"), /Last recorded activity 1m 5s ago/);
  assert.match(run("$('progress-caption').textContent"), /runtime is still running/);
  assert.equal(byClass(run("$('live-timeline')"), 'investigation-update').length, 0);
  run("active.progress.push({kind:'scientific_call',status:'completed',title:'Scientific call: disconnect_target',at:new Date().toISOString()}); updateProgress()");
  assert.doesNotMatch(run("$('progress-caption').textContent"), /Last recorded activity/);
});

test('an activity response begun before submission cannot hide its Stop button', async () => {
  const {ids, run, respond} = harness();
  let resolve;
  respond(route => route === '/activity' ? new Promise(done => { resolve = done; }) : Promise.resolve([]));
  const pending = run('poll()');
  run("activityVersion++; loaded=true; active={conversation_id:'new',status:'queued'}");
  resolve({active:null});
  await pending;
  assert.equal(run('active.conversation_id'), 'new');
  assert.equal(ids.cancel.hidden, false);
});

test('startup discovers the running chat and exposes Stop', async () => {
  const {ids, run, respond} = harness();
  respond(async route => {
    if (route === '/config') return {token:'token',runtime:'test',model:'test'};
    if (route === '/activity') return {active:{conversation_id:'working',status:'running',progress:[]}};
    if (route === '/conversations') return [{id:'working',title:'My investigation'}];
    if (route === '/conversations/working') return {id:'working',title:'My investigation',turns:[{id:'turn',question:'Question',status:'running'}]};
    throw new Error('Unexpected route: ' + route);
  });
  await run('initialize()');
  assert.equal(run('identity'), 'working');
  assert.equal(ids.cancel.hidden, false);
  assert.equal(ids.error.textContent, '');
  assert.equal(ids.saved.children[0].attributes['aria-current'], 'page');
});

test('sending a question opens its chat and replaces Send with live progress and Stop', async () => {
  const {ids, run, respond} = harness();
  run('loaded=true'); ids.question.value = 'Analyze CCO';
  let submitted;
  respond(async (route, options) => {
    if (route === '/turns') {
      submitted = JSON.parse(options.body);
      return {conversation_id:'created',turn_id:'turn'};
    }
    if (route === '/conversations') return [{id:'created',title:'Analyze CCO'}];
    if (route === '/activity') return {active:{conversation_id:'created',status:'running',progress:[]}};
    if (route === '/conversations/created') return {title:'Analyze CCO',turns:[{id:'turn',question:'Analyze CCO',status:'running'}]};
    throw new Error('Unexpected route: ' + route);
  });
  await ids.form.onsubmit({preventDefault() {}});
  assert.deepEqual(submitted, {question:'Analyze CCO',conversation_id:null});
  assert.equal(ids.question.value, '');
  assert.equal(run('identity'), 'created');
  assert.equal(ids.cancel.hidden, false);
  assert.equal(run("$('progress-label').textContent"), 'Investigating…');
});

test('a busy response from another tab reveals that investigation and keeps the draft', async () => {
  const {ids, run, respond} = harness();
  run('loaded=true'); ids.question.value = 'Keep this draft';
  respond(async route => {
    if (route === '/turns') throw new Error('An investigation is running; wait or cancel it first');
    if (route === '/activity') return {active:{conversation_id:'other-tab',status:'running'}};
    if (route === '/conversations') return [{id:'other-tab',title:'Already running'}];
    throw new Error('Unexpected route: ' + route);
  });
  await ids.form.onsubmit({preventDefault() {}});
  assert.equal(ids.question.value, 'Keep this draft');
  assert.equal(ids['global-activity'].hidden, false);
  assert.equal(ids.cancel.hidden, false);
  assert.equal(ids.send.disabled, true);
  assert.match(ids.error.textContent, /investigation is running/);
});
