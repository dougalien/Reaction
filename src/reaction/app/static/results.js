"use strict";

const $ = (id) => document.getElementById(id);
const plot = new LinePlot($("plot"));
let snap = null, result = null, loadingId = null, tab = "species";
let other = null;                            // the run compared against, if any
const W = new URLSearchParams(location.search).get("w") === "2" ? "2" : "1";
const cache = new Map();
function win() { return W === "1" ? { current: snap.current, view: snap.view } : snap.windows["2"]; }
function setView(patch) { return api("/api/view", { patch, window: W }); }
async function load(id) {
  if (!id) return null;
  if (!cache.has(id)) {
    const r = await fetch("/api/run/" + id);
    if (!r.ok) return null;
    cache.set(id, await r.json());
    if (cache.size > 8) cache.delete(cache.keys().next().value);
  }
  return cache.get(id);
}
const MODE_LABEL = { static: "Static", flow_through: "Flow-through" };
const DASH = { static: "", flow_through: "7 5" };
const UNIT_SCALE = { mol: 1, mmol: 1e3, umol: 1e6, mg: 1e3 };
const UNIT_TEXT = { mol: "mol", mmol: "mmol", umol: "µmol", mg: "mg" };
const TABS = [["species", "Species"], ["elements", "Elements"], ["minerals", "Mineral amounts"],
  ["saturation", "Saturation index"], ["property", "Properties"]];
const PROPS = { pH: "pH", ionic_strength: "Ionic strength", water_kg: "Solvent mass", log_aw: "log water activity", log_fO2: "log fO2", Eh: "Eh (V)" };

function pill(el, text, kind) { el.textContent = text; el.className = "pill " + (kind || ""); }

async function onState(s) {
  snap = s;
  const w = win();
  if (w.current && (!result || result.id !== w.current) && loadingId !== w.current) {
    loadingId = w.current;
    try { result = (await load(w.current)) || result; } finally { loadingId = null; }
  }
  if (!w.current) result = null;
  other = w.view.compare ? await load(w.view.compare) : null;
  const run = $("pill-run");
  run.style.display = s.run.state === "running" ? "" : "none";
  pill(run, `calculating ${Math.round(s.run.progress * 100)}%`, "warn");
  draw();
}

connect(onState, (up) => pill($("pill-link"), up ? "connected" : "reconnecting", up ? "ok" : "bad"));

function labelOf(item) {
  const cat = result.catalog;
  if (item.group === "species") return cat.species[item.key].label;
  if (item.group === "elements") return item.key + " (total dissolved)";
  if (item.group === "minerals") return cat.minerals[item.key].label;
  if (item.group === "saturation") return "SI " + cat.minerals[item.key].label;
  return PROPS[item.key];
}

function mwOf(res, item) {
  const c = res.catalog[item.group === "minerals" ? "minerals" : item.group];
  return c && c[item.key] ? c[item.key].mw : null;
}

function valuesOf(item, mode, view, res = result) {
  const run = res.modes[mode];
  if (!run) return null;
  let arr;
  if (item.group === "property") arr = run.properties[item.key];
  else arr = run[item.group][item.key];
  if (!arr) return null;
  if (item.group === "species" || item.group === "elements" || item.group === "minerals") {
    let k = UNIT_SCALE[view.unit];
    if (view.unit === "mg") {
      const mw = mwOf(res, item);
      if (!mw) return null;
      k *= mw;
    }
    return arr.map((v) => (v === null ? null : v * k));
  }
  return arr;
}

function yLabel(items, view) {
  if (!items.length) return "";
  const g = items[0].group;
  if (g === "species" || g === "elements") return `${UNIT_TEXT[view.unit]} per kg of water`;
  if (g === "minerals") return `${UNIT_TEXT[view.unit]} per kg of original water`;
  if (g === "saturation") return "log Q/K (0 = at saturation)";
  return PROPS[items[0].key];
}

function currentSeries() {
  const view = win().view;
  const modes = view.modes.filter((m) => result.modes[m]);
  const series = [];
  view.items.forEach((item, i) => {
    modes.forEach((m) => {
      const y = valuesOf(item, m, view);
      if (!y) return;
      series.push({
        label: labelOf(item) + (modes.length > 1 ? ` (${MODE_LABEL[m].toLowerCase()})` : ""),
        color: `var(--s${(i % 8) + 1})`, dash: DASH[m], y,
        csv: `${labelOf(item)}${modes.length > 1 ? " " + MODE_LABEL[m].toLowerCase() : ""} [${yLabel([item], view)}]`,
      });
    });
    const ox = other && other.axes && other.axes[view.x];
    if (!ox) return;
    // the compared run, in the modes shown here if it has them, otherwise in its own
    const shared = modes.filter((m) => other.modes[m]);
    const omodes = shared.length ? shared : Object.keys(other.modes);
    omodes.forEach((m) => {
      const oy = valuesOf(item, m, view, other);
      if (!oy) return;
      series.push({
        label: `${labelOf(item)} (compared run${omodes.length > 1 || !shared.length ? ", " + MODE_LABEL[m].toLowerCase() : ""})`,
        color: `color-mix(in srgb, var(--s${(i % 8) + 1}) 45%, var(--faint))`, dash: m === "static" ? "2 4" : "9 3 2 3",
        y: oy, x: ox.values, compare: true,
      });
    });
  });
  return series;
}

function drawPickers(view) {
  const runs = snap.runs || [];
  const opt = (r) => h("option", { value: r.id }, `${r.created.slice(5, 16).replace("T", " ")} · ${r.summary}`);
  $("run-pick").replaceChildren(...runs.map(opt));
  $("run-pick").value = win().current || "";
  $("cmp-pick").replaceChildren(h("option", { value: "" }, "Nothing"),
    ...runs.filter((r) => r.id !== win().current).map(opt));
  $("cmp-pick").value = view.compare || "";
  const ox = other && other.axes && other.axes[view.x];
  $("cmp-note").textContent = other && !ox ? `That run has no ${result.axes[view.x].label.toLowerCase()} axis, so it cannot be overlaid here.` : "";
}
$("run-pick").addEventListener("change", (ev) => api("/api/current", { id: ev.target.value, window: W }));
$("cmp-pick").addEventListener("change", (ev) => setView({ compare: ev.target.value || null }));
if (W === "2") { document.title = "Reaction: results 2"; $("open-second").style.display = "none"; }
$("open-second").addEventListener("click", () => window.open("/results?w=2", "reaction-results-2", "popup=no"));

function draw() {
  const has = !!(snap && result);
  $("empty").style.display = has ? "none" : "";
  $("content").style.display = has ? "flex" : "none";
  if (!has) return;
  const view = win().view;
  drawPickers(view);
  $("title").textContent = result.title || result.scenario.title;
  $("meta").textContent = `${result.summary} · ${result.created.replace("T", " ")} · ${result.elapsed_s} s`;
  $("story").replaceChildren(...(result.story || []).map((p) => h("p", {}, p)));
  $("notes").replaceChildren(...result.notes.map((n) => h("li", {}, n)));
  document.title = `Reaction${W === "2" ? " (2)" : ""}: ` + (result.title || result.scenario.title);

  // horizontal axis: only a quantity that changes along the path, and a choice if both do
  const varies = result.varies || ["temperature_c"];
  const xkey = varies.includes(view.x) ? view.x : result.x.key;
  const xax = result.axes ? result.axes[xkey] : result.x;
  const xb = $("xaxis");
  xb.style.display = varies.length > 1 ? "" : "none";
  xb.replaceChildren(...varies.map((k) => h("button", {
    class: k === xkey ? "on" : "", onclick: () => setView({ x: k }) }, result.axes[k].label)));

  const modeIds = Object.keys(result.modes);
  const m = $("modes");
  m.style.display = modeIds.length > 1 ? "" : "none";
  m.replaceChildren(...modeIds.map((id) => h("button", {
    class: view.modes.includes(id) ? "on" : "",
    onclick: () => {
      const next = view.modes.includes(id) ? view.modes.filter((x) => x !== id) : [...view.modes, id];
      if (next.length) setView({ modes: next });
    } }, MODE_LABEL[id])));

  const items = view.items;
  const g = items.length ? items[0].group : "";
  const logOk = g === "species" || g === "elements" || g === "minerals" || (g === "property" && items[0].key === "ionic_strength");
  const useLog = view.log && logOk;
  document.querySelectorAll("#axis button").forEach((b) => {
    b.classList.toggle("on", (b.dataset.log === "1") === useLog);
    b.onclick = () => setView({ log: b.dataset.log === "1" });
  });
  $("axisnote").textContent = view.log && !logOk && items.length ? "Log axis is not used for this quantity" : "";
  const unitsOn = g === "species" || g === "elements" || g === "minerals";
  $("units").style.display = unitsOn ? "" : "none";
  document.querySelectorAll("#units button").forEach((b) => {
    b.classList.toggle("on", b.dataset.unit === view.unit);
    b.onclick = () => setView({ unit: b.dataset.unit });
  });

  const series = currentSeries();
  $("legend").replaceChildren(...series.map((s) => h("span", {},
    h("span", { class: "sw", style: s.dash
      ? `border-color:color-mix(in srgb, ${s.color} 55%, var(--ink));border-top-style:dashed;border-top-width:2px`
      : `border-color:${s.color}` }), s.label)));
  const rel = result.reliability ? result.reliability.levels : [];
  const anyRisk = rel.some((l) => l !== "good");
  const other = varies.length > 1 ? result.axes[varies.find((k) => k !== xkey)] : null;
  $("reliab").style.display = anyRisk ? "" : "none";
  plot.set(xax.values, series, {
    log: useLog, yLabel: yLabel(items, view), xLabel: `${xax.label} (${xax.unit})`,
    xName: xax.label, xUnit: " " + xax.unit, refLines: g === "saturation" ? [0] : null,
    title: result.summary, strip: anyRisk ? rel : null,
    sub: other ? (i) => `${other.label} ${fmtVal(other.values[i])} ${other.unit}` : null,
  });
  drawPicker();
}

function candidates() {
  const run = result.modes[Object.keys(result.modes)[0]];
  const peak = (arr) => Math.max(0, ...arr.filter((v) => v !== null));
  if (tab === "species") {
    return Object.keys(run.species).sort((a, b) => peak(run.species[b]) - peak(run.species[a]))
      .map((k) => ({ group: "species", key: k, label: result.catalog.species[k].label }));
  }
  if (tab === "elements") return Object.keys(run.elements).map((k) => ({ group: "elements", key: k, label: k }));
  if (tab === "minerals") {
    return Object.keys(run.minerals).map((k) => ({ group: "minerals", key: k, label: result.catalog.minerals[k].label }));
  }
  if (tab === "saturation") {
    return Object.keys(run.saturation).map((k) => ({ group: "saturation", key: k,
      label: result.catalog.minerals[k].label + (result.catalog.minerals[k].blocked ? " (blocked)" : "") }));
  }
  return Object.entries(PROPS).filter(([k]) => k in run.properties).map(([k, label]) => ({ group: "property", key: k, label }));
}

function drawPicker() {
  $("tabs").replaceChildren(...TABS.map(([id, label]) => h("button", {
    class: tab === id ? "on" : "", onclick: () => { tab = id; $("search").value = ""; drawPicker(); } }, label)));
  const q = $("search").value.trim().toLowerCase();
  const chosen = win().view.items;
  const list = candidates().filter((c) => !q || c.label.toLowerCase().includes(q) || c.key.toLowerCase().includes(q));
  $("picker").replaceChildren(...list.map((c) => {
    const on = chosen.some((i) => i.group === c.group && i.key === c.key);
    return h("span", { class: "chip" + (on ? " on" : ""), role: "button", tabindex: "0", "aria-pressed": on ? "true" : "false",
      onclick: () => setView({ toggle: { group: c.group, key: c.key } }),
      onkeydown: (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); e.target.click(); } } }, c.label);
  }));
  if (!list.length) $("picker").replaceChildren(h("span", { class: "hint" }, "Nothing matches."));
}
$("search").addEventListener("input", () => snap && result && drawPicker());

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = h("a", { href: url, download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
$("dl-csv").addEventListener("click", () => {
  const series = currentSeries().filter((s) => !s.compare);
  const q = (s) => `"${s.replace(/"/g, '""')}"`;
  const cols = ["reacted", "water_rock", "temperature_c", "pressure_bar"].filter((k) => result.axes[k]).map((k) => result.axes[k]);
  const rows = [[...cols.map((a) => q(`${a.label} (${a.unit})`)), ...series.map((s) => q(s.csv))].join(",")];
  cols[0].values.forEach((_, i) => rows.push([...cols.map((a) => (a.values[i] === null ? "" : a.values[i])),
    ...series.map((s) => (s.y[i] === null ? "" : s.y[i]))].join(",")));
  download(`${result.id}.csv`, rows.join("\n"), "text/csv");
});
$("dl-svg").addEventListener("click", () => download(`${result.id}.svg`, plot.svgText(), "image/svg+xml"));
