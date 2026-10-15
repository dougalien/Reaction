"use strict";

const $ = (id) => document.getElementById(id);
let snap = null, waiting = false, chatCount = 0;

const FLUIDS = [["seawater", "Seawater"]];
FLUIDS.forEach(([k, l]) => $("f-fluid").append(h("option", { value: k }, l)));

function pill(el, text, kind) { el.textContent = text; el.className = "pill " + (kind || ""); }

function render(s) {
  snap = s;
  const c = s.card;
  setIfIdle($("f-title"), c.title);
  $("f-fluid").value = c.fluid;
  $("f-activity_model").value = c.activity_model;
  for (const k of ["pressure_bar", "t_start_c", "t_end_c", "steps"]) setIfIdle($("f-" + k), c[k]);
  document.querySelectorAll("#modeseg button").forEach((b) => b.classList.toggle("on", b.dataset.mode === c.mode));

  const blocked = $("blocked");
  blocked.replaceChildren(...c.blocked_minerals.map((m) => h("span", { class: "chip on static" }, m,
    h("button", { class: "ghost", style: "padding:0 4px", title: "Allow " + m, "aria-label": "Allow " + m,
      onclick: () => api("/api/card", { changes: { blocked_minerals: c.blocked_minerals.filter((x) => x !== m) } }) }, "×"))));
  if (!c.blocked_minerals.length) blocked.replaceChildren(h("span", { class: "hint" }, "none: every mineral in the database may form"));

  $("problems").replaceChildren(...s.card_problems.map((p) => h("div", { class: "warnbox", style: "margin-bottom:4px" }, p)));

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
for (const k of ["title", "pressure_bar", "t_start_c", "t_end_c", "steps"]) {
  $("f-" + k).addEventListener("change", (ev) => edit(k, ev.target.value));
}
$("f-fluid").addEventListener("change", (ev) => edit("fluid", ev.target.value));
$("f-activity_model").addEventListener("change", (ev) => edit("activity_model", ev.target.value));
document.querySelectorAll("#modeseg button").forEach((b) => b.addEventListener("click", () => edit("mode", b.dataset.mode)));
$("add-blocked").addEventListener("keydown", (ev) => {
  if (ev.key !== "Enter") return;
  const v = ev.target.value.trim().toLowerCase();
  if (v && snap && !snap.card.blocked_minerals.includes(v)) edit("blocked_minerals", [...snap.card.blocked_minerals, v]);
  ev.target.value = "";
});
$("run").addEventListener("click", () => api("/api/run"));
$("reset").addEventListener("click", () => api("/api/reset"));
$("open-results").addEventListener("click", () => {
  window.open("/results", "reaction-results", "popup=no");
});
