// Controller tests with a small DOM double: no browser, network, or model calls.
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '../..');

class Node {
  constructor(tag = 'div') {
    this.tag = tag; this.children = []; this.dataset = {}; this.style = {};
    this.attributes = {}; this.value = ''; this.hidden = false; this.disabled = false;
    this.classList = {toggle() {}}; this.scrollHeight = 400; this.scrollTop = 0; this.clientHeight = 400;
  }
  set textContent(value) { this.text = value; this.children = []; }
  get textContent() { return this.text || ''; }
  append(...nodes) { this.children.push(...nodes); }
  replaceChildren(...nodes) { this.children = [...nodes]; }
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

test('saved answer is primary; scientific details and sources start collapsed', () => {
  const {run} = harness();
  const card = run(`answerCard({id:'one',status:'completed',answer:{answer_markdown:'A result',evidence_refs:['sha256:abc'],uncertainties:['Unverified']},
    answer_presentation:{structures:[{smiles:'CCO',image_url:'data:image/svg+xml;base64,abc'}]}})`);
  assert.equal(card.children[0].textContent, 'A result');
  assert.equal(card.querySelectorAll('details').length, 3);
  assert.equal(card.querySelectorAll('details[open]').length, 0);
  assert.ok(card.descendants().some(node => node.tag === 'img' && node.alt.includes('CCO')));
  assert.ok(card.descendants().some(node => node.tag === 'a' && node.textContent === 'Download SVG'));
});

test('progress uses recorded events, elapsed time, and keeps details open across updates', () => {
  const {ids, run} = harness();
  run("identity='working'; active={conversation_id:'working', status:'running', created_at:new Date(Date.now()-65000).toISOString(),progress:[{kind:'web_search',status:'item.completed',at:new Date().toISOString()}]}; renderConversation({title:'Question',turns:[{id:'t1',question:'Question',status:'running'}]})");
  const detail = ids.messages.querySelectorAll('details')[0]; detail.open = true;
  assert.match(run("$('elapsed').textContent"), /1m 5s/);
  assert.equal(run("$('progress-label').textContent"), 'Investigating…');
  assert.ok(ids.messages.descendants().some(node => node.textContent === 'Web search · finished'));
  run("active.progress.push({kind:'command_execution',status:'in_progress'}); renderConversation({title:'Question',turns:[{id:'t1',question:'Question',status:'running'}]})");
  assert.equal(ids.messages.querySelectorAll('details')[0], detail);
  assert.equal(detail.open, true);
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
