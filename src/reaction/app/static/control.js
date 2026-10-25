"use strict";

const $ = (id) => document.getElementById(id);
let snap = null, waiting = false, chatCount = 0;

let fluidKeys = "";

function pill(el, text, kind) { el.textContent = text; el.className = "pill " + (kind || ""); }

function render(s) {
  snap = s;
  const c = s.card;
  setIfIdle($("f-title"), c.title);
  renderFluid(s);
  $("f-activity_model").value = c.activity_model;
  for (const k of ["p_start_bar", "p_end_bar", "t_start_c", "t_end_c", "steps"]) setIfIdle($("f-" + k), c[k]);
  document.querySelectorAll("#modeseg button").forEach((b) => b.classList.toggle("on", b.dataset.mode === c.mode));

  const reactants = $("reactants");
  reactants.replaceChildren(...c.reactants.map((r, i) => h("span", { class: "chip on static" },
    `${r.amount} ${r.unit} ${r.name}`,
    h("button", { class: "ghost", style: "padding:0 4px", title: "Remove " + r.name, "aria-label": "Remove " + r.name,
      onclick: () => edit("reactants", c.reactants.filter((_, k) => k !== i)) }, "×"))));
  if (!c.reactants.length) reactants.replaceChildren(h("span", { class: "hint" }, "none: only temperature and pressure change"));

  const blocked = $("blocked");
  blocked.replaceChildren(...c.blocked_minerals.map((m) => h("span", { class: "chip on static" }, m,
    h("button", { class: "ghost", style: "padding:0 4px", title: "Allow " + m, "aria-label": "Allow " + m,
      onclick: () => api("/api/card", { changes: { blocked_minerals: c.blocked_minerals.filter((x) => x !== m) } }) }, "×"))));
  if (!c.blocked_minerals.length) blocked.replaceChildren(h("span", { class: "hint" }, "none: every mineral in the database may form"));

  $("problems").replaceChildren(...s.card_problems.map((p) => h("div", { class: "errbox", style: "margin-bottom:4px" }, p)));
  $("advice").replaceChildren(...(s.advisories || []).map((p) => h("div", { class: "warnbox", style: "margin-bottom:4px" }, p)));

  const running = s.run.state === "running";
  $("run").disabled = running;
  $("run").textContent = running ? "Running…" : "Run";
  $("prog").style.display = running ? "block" : "none";
  $("prog").value = s.run.progress;
  $("runmsg").textContent = running ? s.run.message : "";
  const err = $("runerr");
  err.style.display = s.run.state === "error" ? "block" : "none";
  err.textContent = s.run.error;

  const chat = $("chat");
  const atBottom = chat.scrollHeight - chat.scrollTop - chat.clientHeight < 80;
  if (s.chat.length !== chatCount) waiting = waiting && s.chat[s.chat.length - 1].role === "user";
  chatCount = s.chat.length;
  const bubbles = s.chat.map((m) => h("div", { class: "bub " + m.role }, m.text));
  if (waiting) bubbles.push(h("div", { class: "bub assistant thinking" }, "Thinking…"));
  chat.replaceChildren(...bubbles);
  if (atBottom || waiting) chat.scrollTop = chat.scrollHeight;

  $("hist").replaceChildren(...s.runs.map((r) => h("button", {
    class: r.id === s.current ? "on" : "", onclick: () => api("/api/current", { id: r.id }),
    title: r.summary }, h("div", {}, r.summary), h("div", { class: "hint" }, r.created.replace("T", " ")))));
  if (!s.runs.length) $("hist").replaceChildren(h("div", { class: "hint" }, "No runs yet."));

  pill($("pill-llm"), s.llm.ready ? `${s.llm.model} ready` : `${s.llm.model || "language model"} offline: simple rules only`,
    s.llm.ready ? "ok" : "warn");
  const e = s.engine;
  pill($("pill-engine"), e.ready === false ? "data missing" : e.ready ? "data loaded" : "data loads on first run",
    e.ready === false ? "bad" : e.ready ? "ok" : "");
  $("pill-engine").title = e.message;
}

function renderFluid(s) {
  const c = s.card, read = s.fluid_read || { lines: [] }, fluids = s.fluids || [];
  const keys = fluids.map((f) => f.key).join();
  if (keys !== fluidKeys) {
    fluidKeys = keys;
    $("f-fluid").replaceChildren(...fluids.map((f) => h("option", { value: f.key }, f.label)));
    $("f-mix_with").replaceChildren(h("option", { value: "" }, "None"),
      ...fluids.filter((f) => f.key !== "custom").map((f) => h("option", { value: f.key }, f.label)));
  }
  $("f-mix_with").value = c.mix_with;
  for (const k of ["mix_ratio", "mix_t_c"]) {
    $("f-" + k).disabled = !c.mix_with;
    setIfIdle($("f-" + k), c[k] === null ? "" : c[k]);
  }
  $("mix-hint").textContent = !c.mix_with ? "" : c.mix_t_c === null
    ? "Mixed in steadily; the temperature follows the path above."
    : "Mixed in steadily; the temperature follows the mixing (mass-weighted, equal heat capacities), so the end temperature above is not used.";
  $("f-t_end_c").disabled = c.mix_with && c.mix_t_c !== null;
  $("f-fluid").value = c.fluid;
  const custom = c.fluid === "custom";
  const f = fluids.find((x) => x.key === c.fluid);
  $("fluid-desc").textContent = custom ? "" :
    (f ? f.description + " Editing it below starts a custom fluid from it." : "");
  $("custom").style.display = "block";
  {
    setIfIdle($("f-analysis"), read.analysis);
    $("f-units").value = read.units;
    setIfIdle($("f-ph"), read.ph === null ? "" : read.ph);
    setIfIdle($("f-fluid_t_c"), read.measured_t_c);
    const ions = read.lines.filter((l) => l.charge !== 0).map((l) => l.species);
    const opts = [["auto", "Auto"], ["none", "None (use as given)"], ...ions.map((i) => [i, "On " + i])];
    if (!opts.some(([v]) => v === read.charge_balance)) opts.push([read.charge_balance, "On " + read.charge_balance]);
    $("f-charge_balance").replaceChildren(...opts.map(([v, l]) => h("option", { value: v }, l)));
    $("f-charge_balance").value = read.charge_balance;
  }
  const n = read.lines.length;
  $("fluid-sum").textContent = n ? `As read: ${n} solute${n > 1 ? "s" : ""}, pH ${read.ph ?? "not given"} at ${read.measured_t_c} °C, ` +
    `${fmtVal(read.tds_mg_per_kg)} mg/kg dissolved solids` : "Nothing read yet";
  const rows = read.lines.map((l) => h("tr", {}, h("td", {}, l.written), h("td", {}, `${l.value} ${l.unit}`),
    h("td", { class: "muted" }, l.molal === null ? "" : fmtVal(l.molal) + " mol/kg water")));
  $("fluid-table").replaceChildren(h("table", { class: "small", style: "margin-top:4px;border-spacing:10px 1px" }, ...rows),
    ...read.notes.map((t) => h("div", { class: "hint" }, t)));
}

connect(render, (up) => pill($("pill-link"), up ? "connected" : "reconnecting", up ? "ok" : "bad"));

$("compose").addEventListener("submit", (ev) => {
  ev.preventDefault();
  const text = $("msg").value.trim();
  if (!text) return;
  $("msg").value = "";
  waiting = true;
  render({ ...snap, chat: [...snap.chat, { role: "user", text, t: "" }] });
  api("/api/chat", { text });
});
$("msg").addEventListener("keydown", (ev) => {
  if (ev.key === "Enter" && !ev.shiftKey) { ev.preventDefault(); $("compose").requestSubmit(); }
});

function edit(field, value) { api("/api/card", { changes: { [field]: value } }); }
for (const k of ["title", "p_start_bar", "p_end_bar", "t_start_c", "t_end_c", "steps"]) {
  $("f-" + k).addEventListener("change", (ev) => edit(k, ev.target.value));
}
$("f-fluid").addEventListener("change", (ev) => edit("fluid", ev.target.value));
for (const k of ["analysis", "units", "ph", "fluid_t_c", "charge_balance", "mix_with", "mix_ratio", "mix_t_c"]) {
  $("f-" + k).addEventListener("change", (ev) => edit(k, ev.target.value));
}
$("f-activity_model").addEventListener("change", (ev) => edit("activity_model", ev.target.value));
document.querySelectorAll("#modeseg button").forEach((b) => b.addEventListener("click", () => edit("mode", b.dataset.mode)));
$("add-blocked").addEventListener("keydown", (ev) => {
  if (ev.key !== "Enter") return;
  const v = ev.target.value.trim().toLowerCase();
  if (v && snap && !snap.card.blocked_minerals.includes(v)) edit("blocked_minerals", [...snap.card.blocked_minerals, v]);
  ev.target.value = "";
});
$("add-reactant").addEventListener("keydown", async (ev) => {
  if (ev.key !== "Enter") return;
  const v = ev.target.value.trim();
  if (!v || !snap) return;
  const reply = await api("/api/card", { changes: { reactants: [...snap.card.reactants.map((r) => `${r.amount} ${r.unit} ${r.name}`), v].join("; ") } });
  if (!(reply.problems || []).length) ev.target.value = "";
  else $("add-reactant").setCustomValidity(reply.problems.join("; ")), $("add-reactant").reportValidity();
});
$("add-reactant").addEventListener("input", (ev) => ev.target.setCustomValidity(""));
$("run").addEventListener("click", () => api("/api/run"));
$("reset").addEventListener("click", () => api("/api/reset"));
$("open-results").addEventListener("click", () => {
  window.open("/results", "reaction-results", "popup=no");
});
