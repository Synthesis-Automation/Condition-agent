// Controller tests with a small DOM double: no browser, network, or model calls.
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '../..');

test('prepared literature cards distinguish indexed graphs, acquisition and conflicting evidence', () => {
  const {run, context} = harness();
  context.prepared = {
    title:'Source reaction', source_id:'paper', locator:'Example 8', structure_origin:'indexed_record',
    image_url:'data:image/svg+xml;base64,source', reaction_smiles:'CCBr.N>>CCN',
    structure_evidence:'Compound 8 gave compound 9.', conditions:[], limitations:['Unresolved source discrepancy: DCM versus DCE.'],
    source_provenance:{acquisition:'agent_supplied_excerpt',extraction_status:'agent_supplied',limitations:[]},
    reactants:[{name:'Bromoethane',compound_id:'8',material_form:'Neutral parent; concentration unknown',
      evidence_artifact_url:'/saved/reactant',graph_status:'graph_checked_assignment_unverified',formula:'C2H5Br'}],
    products:[{name:'Product',compound_id:'9',material_form:'Salt form unknown',graph_status:'conflicting',
      evidence_artifact_url:'/saved/product',formula:'C2H7N'}],
    preparation_artifact_url:'/saved/preparation'
  };
  const before=JSON.stringify(context.prepared);
  const card=run("literatureReactionCard(prepared, new Map([['paper',{title:'Paper',url:'https://example.org/paper'}]]), 'lit')");
  const texts=card.descendants().map(node=>node.textContent);
  assert.ok(texts.includes('Indexed source structures'));
  assert.ok(texts.some(text=>text.includes('acquisition and transcription are not independently verified')));
  assert.ok(texts.some(text=>text.includes('Source and graph evidence conflict')));
  assert.ok(texts.includes('Unresolved source discrepancy: DCM versus DCE.'));
  assert.ok(card.querySelectorAll('a').some(node=>node.href==='/saved/product'));
  assert.ok(card.querySelectorAll('a').some(node=>node.href==='/saved/preparation'));
  assert.equal(JSON.stringify(context.prepared),before);
});

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

test('workspace mode is selectable for new chats and locked to a reopened conversation', async () => {
  const {ids, run, respond} = harness();
  run('loaded=true; updateControls()');
  assert.equal(ids['workspace-mode'].value, 'normal');
  assert.equal(ids['workspace-mode'].disabled, false);
  ids['workspace-mode'].value = 'pure_agent';
  ids['workspace-mode'].onchange();
  run("identity='saved'; renderConversation({title:'Saved',workspace_mode:'tools_only',turns:[]})");
  assert.equal(ids['workspace-mode'].value, 'tools_only');
  assert.equal(ids['workspace-mode'].disabled, true);
  assert.match(ids['workspace-mode'].title, /new chat to change mode/);
  run("renderConversation({title:'Legacy',turns:[]})");
  assert.equal(ids['workspace-mode'].value, 'normal');
  respond(async () => []);
  run('newChat()');
  assert.equal(ids['workspace-mode'].value, 'pure_agent');
  assert.equal(ids['workspace-mode'].disabled, false);
});

for (const mode of ['tools_only', 'tools_formatting', 'tools_guidance']) {
test(`submitting preserves and locks the selected mode: ${mode}`, async () => {
  const {ids, run, respond} = harness();
  run(`loaded=true; newMode='${mode}'; updateControls()`);
  ids.question.value = 'Compare the chemistry';
  let sent;
  respond(async (route, options) => {
    if (route === '/turns') {
      sent = JSON.parse(options.body);
      assert.equal(ids['workspace-mode'].disabled, true);
      return {conversation_id:'chat',turn_id:'turn'};
    }
    if (route === '/conversations/chat') return {id:'chat',title:'Comparison',workspace_mode:mode,turns:[{
      id:'turn',question:'Compare the chemistry',status:'completed',answer:{schema_version:'agent_text.v1',answer_markdown:'Free text.'},
    }]};
    if (route === '/conversations') return [];
    if (route === '/activity') return {active:null};
    throw new Error('Unexpected route: ' + route);
  });
  await ids.form.onsubmit({preventDefault() {}});
  assert.deepEqual(sent, {question:'Compare the chemistry',conversation_id:null,mode});
  assert.equal(ids['workspace-mode'].value, mode);
  assert.equal(ids['workspace-mode'].disabled, true);
  assert.ok(ids.messages.descendants().some(node => node.textContent === 'Free text.'));
});
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

test('step assessments keep structural support separate from recipe coverage', () => {
  const {run, context} = harness();
  context.stepEvidence = {assessment_evidence:{status:'recorded',structural_assessments:[{
    status:'precedent_supported',artifact_url:'/saved/structure',
    gates:[{gate_id:'condition_support',status:'not_run',summary:'No condition evaluation requested.'}],
  }],recipe_assessments:[{status:'no_known_conflict',recipe_id:'recipe-one',artifact_url:'/saved/recipe',
    coverage:{capability_status:'not_covered',condition_identity_status:'resolved',
      evaluated_hard_conflict_rule_ids:['a','b'],evaluated_soft_penalty_rule_ids:['c'],
      evaluated_regime_requirement_ids:[],limitations:['Outcome remains untested.']}}]}};
  let view = run("stepAssessmentEvidence(stepEvidence, 's1')");
  let texts = view.descendants().map(node => node.textContent);
  assert.ok(texts.includes('Structural assessment: precedent supported.'));
  assert.ok(texts.includes('No applicable reaction capability requirement was checked.'));
  assert.ok(texts.includes('Rules evaluated: 2 conflict checks, 1 penalty checks, 0 regime requirements.'));
  assert.ok(texts.includes('condition support: not run — No condition evaluation requested.'));
  assert.ok(texts.includes('Experimental feasibility is not established by these checks.'));
  assert.ok(view.querySelectorAll('a').some(node => node.href === '/saved/recipe'));
  assert.equal(view.querySelectorAll('details[open]').length, 0);
  delete context.stepEvidence.assessment_evidence.recipe_assessments[0].coverage;
  view = run("stepAssessmentEvidence(stepEvidence, 's1')");
  assert.ok(view.descendants().some(node => node.textContent === 'Reaction capability coverage was not recorded.'));
  context.stepEvidence.assessment_evidence = {status:'evidence_unavailable'};
  view = run("stepAssessmentEvidence(stepEvidence, 's1')");
  assert.ok(view.descendants().some(node => node.textContent.includes('support remains unresolved')));
  assert.ok(!view.descendants().some(node => node.textContent.includes('precedent supported')));
});

test('supporting reactions retain source identity, separate observations and evidence gaps', () => {
  const {run, context} = harness();
  const step = {id:'s1',title:'Proposed step',basis:'proposed',reactant_ids:['a'],product_ids:['b'],
    after_step_ids:[],conditions:[],yield_info:null,source_ids:[],limitations:[]};
  context.supportView = {sources:[],molecules:[{id:'a',name:'A'},{id:'b',name:'B'}],routes:[],steps:[step]};
  let card = run('showStructured(supportView)');
  assert.ok(card.descendants().some(node => node.textContent === 'No inspected experimental support supplied for this step.'));
  step.supporting_evidence = [{artifact_ref:'sha256:fixture',artifact_url:'/saved/inspection',status:'precedents_available',
    scope:'selected_template_top_20',saved_match_count:4,distinct_references_on_page:1,page:{next_offset:1},
    precedents:[{match_id:'match1',reaction_id:'reaction-1',reference_title:'Reported patent example',reference_url:'https://example.org/patent',
      reaction_smiles:'CCO>>CC=O',image_url:'data:image/svg+xml;base64,abc',same_recorded_product:false,same_recorded_precursors:false,
      product_similarity:0.4,precursor_similarity:0.3,limitations:['Transfer remains uncertain.'],
      observations:[{observation_id:'obs1',yield_pct:40,resolved_recipe:{solvents:[{canonical_name:'Ethanol'}],warnings:['CONDITION_IDENTITY_UNCERTAINTY']}},{observation_id:'obs2',yield_pct:null}],
      procedures:[{observation_id:'obs1',procedure_text:'Literal source <script>text</script>'}]}]}];
  card = run('showStructured(supportView)');
  assert.equal(card.querySelectorAll('img').length,1);
  assert.ok(card.descendants().some(node => node.textContent === 'Precedent support'));
  assert.equal(byClass(card, 'precedent-card')[0].parentNode.className, 'step-precedents');
  assert.ok(card.querySelectorAll('a').some(node => node.href === 'https://example.org/patent'));
  assert.ok(card.descendants().some(node => node.textContent === 'Reported yield: 40%'));
  assert.ok(card.descendants().some(node => node.textContent === 'Conditions and yield not recorded for this experiment.'));
  const otherExperiment = card.querySelectorAll('summary').find(node => node.textContent === 'Experiment 2');
  assert.ok(!otherExperiment.parentNode.open, 'additional experiments remain separate and collapsed');
  assert.ok(card.querySelectorAll('code').some(node => node.textContent === 'CCO>>CC=O'));
  const conditions = card.descendants().find(node => node.className === 'precedent-conditions');
  assert.equal(conditions.parentNode.className, 'precedent-card');
  assert.ok(conditions.descendants().some(node => node.textContent === 'Reagents & solvents: Ethanol'));
  const cautions = card.querySelectorAll('summary').find(node => node.textContent === 'Match details & cautions');
  assert.ok(!cautions.parentNode.open);
  assert.ok(cautions.parentNode.descendants().some(node=>node.textContent==='CONDITION_IDENTITY_UNCERTAINTY'));
  assert.ok(card.descendants().some(node => node.textContent === 'Literal source <script>text</script>'));
  assert.ok(!card.querySelectorAll('a').some(node => node.download));
  step.supporting_evidence = [{status:'no_precedents_retrieved',precedents:[],scope:'saved_assessment_matches',saved_match_count:0,distinct_references_on_page:0}];
  card = run('showStructured(supportView)');
  assert.ok(card.descendants().some(node => node.textContent === 'No supporting experiment found in this recorded search.'));
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
  assert.equal(images[0].style.width, '735px', 'keep source and target schemes readable at the same scale');
  assert.match(images[0].alt, /Reactant → Product/);
  assert.equal(images[0].src, context.fixture.structured_presentation.steps[0].image_url);
  assert.equal(card.querySelectorAll('a').some(node => node.download), false,
    'reaction schemes should not show a download button');
  const stepDetails = card.querySelectorAll('details').find(node => node.dataset.key === 'scheme:science:s1');
  assert.equal(stepDetails.firstElementChild.textContent, 'Step cautions & assessment details');
  assert.ok(stepDetails.descendants().some(node => node.textContent === 'Feasibility unverified'),
    'material cautions remain accessible beside the step');
  assert.ok(byClass(view, 'step-conditions')[0].descendants().some(node => node.textContent === 'Unknown oxidant'));
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
  assert.ok(alternatives[0].descendants().some(node => node.textContent === 'Supply unconfirmed'));
  assert.equal(Boolean(footer.open), false);
  assert.equal(routedCard.children.filter(node => node.tag === 'details').length, 2);
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
  assert.ok(list.some(node => node.textContent.endsWith(': Stock unknown')));
  assert.ok(list.some(node => node.textContent === 'Route A: Specific precursor needed'));
  assert.ok(byClass(card, 'answer-uncertainties')[0].descendants().some(node => node.textContent === '<script>plain text</script>'));
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
  assert.deepEqual(submitted, {question:'Analyze CCO',conversation_id:null,mode:'normal'});
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


test('linear routes show source SVGs and actual conditions with diagnostics collapsed', () => {
  const {run, context} = harness();
  const attribution = {basis:'proposed',source_ids:[],limitations:['Target feasibility unresolved']};
  const step = {id:'s1',title:'First step',...attribution,reactant_ids:['a'],product_ids:['b'],after_step_ids:[],
    image_url:'data:image/svg+xml;base64,target-step',conditions:[{text:'Target recipe',...attribution}],
    yield_info:{text:'Yield not established',basis:'unknown',source_ids:[],limitations:[]},
    rationale:{text:'Reason for this step',...attribution},
    supporting_evidence:[{status:'precedents_available',saved_match_count:1,precedents:[{
      reaction_id:'p1',reference_title:'Source experiment',image_url:'data:image/svg+xml;base64,source',
      observations:[{yield_pct:72}],compatibility:{status:'unknown',analysis_warnings:['Transfer remains uncertain']},
    }]}]};
  context.routeView = {sources:[],molecules:[{id:'a',name:'A'},{id:'b',name:'B'},{id:'c',name:'C'}],
    steps:[step,{...step,id:'s2',title:'Second step',reactant_ids:['b'],product_ids:['c'],after_step_ids:['s1'],supporting_evidence:[]}],
    routes:[{id:'r1',title:'Complete route',step_ids:['s1','s2'],limitations:['Route is proposed'],unreached_target_ids:[],
      image_url:'data:image/svg+xml;base64,route',scheme_width:1100,drawing_status:'linear_scheme'}]};
  const original = JSON.stringify(context.routeView);
  const card = run('showStructured(routeView)');
  const images = card.querySelectorAll('img');
  assert.equal(images.filter(node => node.src.endsWith('route')).length, 1);
  assert.equal(images.filter(node => node.src.endsWith('target-step')).length, 0);
  const sourceFigure = images.find(node => node.src.endsWith('source')).parentNode.parentNode;
  assert.equal(sourceFigure.parentNode.className, 'precedent-card');
  assert.ok(card.querySelectorAll('figcaption').some(node => node.textContent.includes('left to right')));
  const inDetails = node => {
    for (let parent = node.parentNode; parent; parent = parent.parentNode) if (parent.tag === 'details') return true;
    return false;
  };
  assert.ok(!inDetails(sourceFigure), 'source reaction SVG must be visible');
  for (const text of ['Reason for this step','Target recipe','No inspected experimental support supplied for this step.',
    'Reported yield: 72%']) {
    const nodes = card.descendants().filter(node => node.textContent === text);
    assert.ok(nodes.some(node => !inDetails(node)), text + ' must remain visible');
  }
  for (const text of ['Transfer remains uncertain','Target feasibility unresolved']) {
    assert.ok(card.descendants().filter(node => node.textContent === text).every(inDetails), text + ' is retained in cautions');
  }
  assert.ok(!card.descendants().some(node => node.textContent.startsWith('Conditions: Proposed')));
  assert.equal(JSON.stringify(context.routeView), original);
  context.routeView.routes[0].drawing_status = 'dependency_overview';
  const fallback = run('showStructured(routeView)');
  assert.equal(fallback.querySelectorAll('img').filter(node => node.src.endsWith('target-step')).length, 2);
});

test('condition choices share a scheme and expose two precedents with attributed rationale', () => {
  const {run, context} = harness();
  const claim = {basis:'proposed',source_ids:['paper'],limitations:[]};
  const first = {id:'s1',title:'Preferred recipe',...claim,reactant_ids:['a'],product_ids:['b'],
    after_step_ids:[],reaction_smiles:'CCO>>CC=O',image_url:'data:image/svg+xml;base64,target',scheme_width:980,
    conditions:[{text:'Recorded recipe adaptation',...claim}],yield_info:null,
    rationale:{text:'Preserve the sensitive group <script>literal</script>',...claim},
    supporting_evidence:[{status:'precedents_available',saved_match_count:2,precedents:[
      {reaction_id:'r1',reference_title:'Experiment one',support_kind:'condition_observation',reaction_smiles:'CCO>>CC=O',
        image_url:'data:image/svg+xml;base64,source',observations:[{observation_id:'o1',yield_pct:42}],
        structural_comparison:{environments:{shared:[],query_only:['a'],precedent_only:['b']}},
        compatibility:{status:'unknown',analysis_warnings:['Structural coverage incomplete'],unresolved_requirements:['Catalyst unresolved']}},
      {reaction_id:'r2',reference_title:'Experiment two',reaction_smiles:'CCO>>CC=O',image_url:'data:image/svg+xml;base64,second-source',observations:[]},
      {reaction_id:'r3',reference_title:'Experiment three',reaction_smiles:'CCO>>CC=O',image_url:'data:image/svg+xml;base64,third-source',observations:[]},
    ]}]};
  context.choices={sources:[{id:'paper',kind:'external_source',title:'Captured Example 2',url:'https://example.org/paper',locator:'Example 2'}],
    molecules:[{id:'a',name:'Alcohol'},{id:'b',name:'Aldehyde'}],routes:[],
    steps:[first,{...first,id:'s2',title:'Alternative recipe',source_ids:[],conditions:[{text:'Alternative solvent',...claim}],supporting_evidence:[]}]};
  const card=run('showStructured(choices)');
  assert.equal(card.querySelectorAll('img').filter(node=>node.src.endsWith('target')).length,1);
  const alternatives=card.querySelectorAll('summary').find(node=>node.textContent==='Alternative conditions (1)').parentNode;
  assert.ok(!alternatives.open);
  assert.ok(alternatives.descendants().some(node=>node.textContent==='Alternative solvent'));
  assert.ok(!alternatives.descendants().some(node=>node.textContent==='No inspected experimental support supplied for this step.'));
  assert.ok(alternatives.descendants().some(node=>node.textContent==='Source reaction structures were not supplied for these literature citations.'));
  const support=byClass(card,'step-precedents')[0];
  assert.equal(byClass(support,'precedent-card')[0].parentNode,support);
  assert.equal(byClass(support,'precedent-card')[1].parentNode,support);
  assert.ok(!card.querySelectorAll('summary').find(node=>node.textContent==='More precedents (1)').parentNode.open);
  assert.ok(support.descendants().some(node=>node.textContent==='Recorded differences: environments.'));
  assert.ok(support.descendants().some(node=>node.textContent==='Structural coverage incomplete'));
  assert.ok(support.descendants().some(node=>node.textContent==='Catalyst unresolved'));
  const rationale=byClass(card,'step-rationale')[0];
  assert.ok(rationale.descendants().some(node=>node.textContent==='Proposed'));
  assert.ok(rationale.querySelectorAll('a').some(node=>node.href==='https://example.org/paper'));
  assert.equal(card.querySelectorAll('script').length,0);
  for(const code of card.querySelectorAll('code')) {
    let parent=code.parentNode;
    while(parent && parent.tag!=='details') parent=parent.parentNode;
    assert.ok(parent && !parent.open,'source notation stays in details');
  }
  context.choices.steps[1].reaction_smiles='CCCO>>CCC=O';
  assert.equal(run('showStructured(choices)').querySelectorAll('img').filter(node=>node.src.endsWith('target')).length,2,
    'different explicit transformations must not share a target scheme');
});

test('retrosynthesis leads with chemistry, retains cautions, and exposes requests for user input', () => {
  const {run, context} = harness();
  const claim = {basis:'proposed',source_ids:['paper'],limitations:[]};
  context.chemistryTurn = {id:'retro',status:'completed',progress:[{kind:'agent_update',detail:'Research log'}],
    answer:{answer_markdown:'Long written explanation',uncertainties:['Select a solvent'],evidence_refs:[]},
    structured_presentation:{sources:[{id:'paper',kind:'external_source',title:'Paper Example 3',url:'https://example.org/paper',locator:'Example 3'}],
      molecules:[{id:'a',name:'A',limitations:[]},{id:'b',name:'B',limitations:[]}],claims:[],
      routes:[{id:'route',title:'Route',step_ids:['s'],image_url:'route.svg',drawing_status:'linear_scheme',limitations:[],unreached_target_ids:[]}],
      steps:[{id:'s',title:'Substitution',...claim,reactant_ids:['a'],product_ids:['b'],after_step_ids:[],
        conditions:[{text:'EtOH, 80 °C, 2 h',...claim}],reagents:[{text:'Base',...claim}],yield_info:null,
        rationale:{text:'Replace the leaving group',...claim},
        assessment_evidence:{status:'recorded',structural_assessments:[],recipe_assessments:[{
          status:'conflicting',hard_conflicts:['Known incompatible condition'],coverage:{capability_status:'not_covered'},
        }]},
        supporting_evidence:[{status:'precedents_available',precedents:[{
          reaction_id:'ref',reference_title:'Paper Example 3',reference_url:'https://example.org/paper',image_url:'source.svg',observations:[],
        }]}],
      }],
    }};
  const original = JSON.stringify(context.chemistryTurn);
  const inClosedDetails = node => {
    for (let parent=node.parentNode;parent;parent=parent.parentNode) if (parent.tag==='details' && !parent.open) return true;
    return false;
  };
  const card = run('answerCard(chemistryTurn)');
  assert.equal(card.children[0].className, 'scientific-view');
  for (const text of ['Long written explanation','Select a solvent','Research log','No applicable reaction capability requirement was checked.']) {
    const nodes = card.descendants().filter(n=>n.textContent===text);
    assert.ok(nodes.length && nodes.every(inClosedDetails), text + ' stays in details');
  }
  for (const text of ['Replace the leaving group','EtOH, 80 °C, 2 h','Known incompatible condition']) {
    assert.ok(card.descendants().some(n=>n.textContent===text && !inClosedDetails(n)), text + ' must be visible');
  }
  assert.ok(card.querySelectorAll('img').some(n=>n.src==='source.svg' && !inClosedDetails(n)));
  assert.ok(card.querySelectorAll('a').some(n=>n.href==='https://example.org/paper' && !inClosedDetails(n)));
  assert.equal(byClass(card,'step-conditions')[0].querySelectorAll('a').length, 0, 'source links are collected in support, not repeated per condition');
  assert.equal(JSON.stringify(context.chemistryTurn), original);
  context.chemistryTurn.answer.needs_user_input = true;
  const needsInput = run('answerCard(chemistryTurn)');
  assert.ok(needsInput.descendants().some(n=>n.textContent==='Select a solvent' && !inClosedDetails(n)));
  assert.ok(needsInput.children.some(n=>n.textContent==='Long written explanation'));
});
