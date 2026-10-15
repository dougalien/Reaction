"use strict";

const $ = (id) => document.getElementById(id);
const plot = new LinePlot($("plot"));
let snap = null, result = null, loadingId = null, tab = "species";
const MODE_LABEL = { static: "Static", flow_through: "Flow-through" };
const DASH = { static: "", flow_through: "7 5" };
const UNIT_SCALE = { mol: 1, mmol: 1e3, umol: 1e6 };
const UNIT_TEXT = { mol: "mol", mmol: "mmol", umol: "µmol" };
const TABS = [["species", "Species"], ["elements", "Elements"], ["minerals", "Mineral amounts"],
  ["saturation", "Saturation index"], ["property", "Properties"]];
const PROPS = { pH: "pH", ionic_strength: "Ionic strength", water_kg: "Solvent mass", log_aw: "log water activity" };

function pill(el, text, kind) { el.textContent = text; el.className = "pill " + (kind || ""); }

async function onState(s) {
  snap = s;
  if (s.current && (!result || result.id !== s.current) && loadingId !== s.current) {
    loadingId = s.current;
    try {
      const r = await fetch("/api/run/" + s.current);
      if (r.ok) result = await r.json();
    } finally { loadingId = null; }
  }
  if (!s.current) result = null;
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

function valuesOf(item, mode, view) {
  const run = result.modes[mode];
  let arr;
  if (item.group === "property") arr = run.properties[item.key];
  else arr = run[item.group][item.key];
  if (!arr) return null;
  if (item.group === "species" || item.group === "elements" || item.group === "minerals") {
    const k = UNIT_SCALE[view.unit];
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
  const view = snap.view;
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
  });
  return series;
}

function draw() {
  const has = !!(snap && result);
  $("empty").style.display = has ? "none" : "";
  $("content").style.display = has ? "flex" : "none";
  if (!has) return;
  const view = snap.view;
  $("title").textContent = result.title || result.scenario.title;
  $("meta").textContent = `${result.summary} · ${result.created.replace("T", " ")} · ${result.elapsed_s} s`;
  $("story").replaceChildren(...(result.story || []).map((p) => h("p", {}, p)));
  $("notes").replaceChildren(...result.notes.map((n) => h("li", {}, n)));
  document.title = "Reaction: " + (result.title || result.scenario.title);

  // horizontal axis: only a quantity that changes along the path, and a choice if both do
  const varies = result.varies || ["temperature_c"];
  const xkey = varies.includes(view.x) ? view.x : result.x.key;
  const xax = result.axes ? result.axes[xkey] : result.x;
  const xb = $("xaxis");
  xb.style.display = varies.length > 1 ? "" : "none";
  xb.replaceChildren(...varies.map((k) => h("button", {
    class: k === xkey ? "on" : "", onclick: () => api("/api/view", { patch: { x: k } }) }, result.axes[k].label)));

  const modeIds = Object.keys(result.modes);
  const m = $("modes");
  m.style.display = modeIds.length > 1 ? "" : "none";
  m.replaceChildren(...modeIds.map((id) => h("button", {
    class: view.modes.includes(id) ? "on" : "",
    onclick: () => {
      const next = view.modes.includes(id) ? view.modes.filter((x) => x !== id) : [...view.modes, id];
      if (next.length) api("/api/view", { patch: { modes: next } });
    } }, MODE_LABEL[id])));

  const items = view.items;
  const g = items.length ? items[0].group : "";
  const logOk = g === "species" || g === "elements" || g === "minerals" || (g === "property" && items[0].key === "ionic_strength");
  const useLog = view.log && logOk;
  document.querySelectorAll("#axis button").forEach((b) => {
    b.classList.toggle("on", (b.dataset.log === "1") === useLog);
    b.onclick = () => api("/api/view", { patch: { log: b.dataset.log === "1" } });
  });
  $("axisnote").textContent = view.log && !logOk && items.length ? "Log axis is not used for this quantity" : "";
  const unitsOn = g === "species" || g === "elements" || g === "minerals";
  $("units").style.display = unitsOn ? "" : "none";
  document.querySelectorAll("#units button").forEach((b) => {
    b.classList.toggle("on", b.dataset.unit === view.unit);
    b.onclick = () => api("/api/view", { patch: { unit: b.dataset.unit } });
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
  return Object.entries(PROPS).map(([k, label]) => ({ group: "property", key: k, label }));
}

function drawPicker() {
  $("tabs").replaceChildren(...TABS.map(([id, label]) => h("button", {
    class: tab === id ? "on" : "", onclick: () => { tab = id; $("search").value = ""; drawPicker(); } }, label)));
  const q = $("search").value.trim().toLowerCase();
  const chosen = snap.view.items;
  const list = candidates().filter((c) => !q || c.label.toLowerCase().includes(q) || c.key.toLowerCase().includes(q));
  $("picker").replaceChildren(...list.map((c) => {
    const on = chosen.some((i) => i.group === c.group && i.key === c.key);
    return h("span", { class: "chip" + (on ? " on" : ""), role: "button", tabindex: "0", "aria-pressed": on ? "true" : "false",
      onclick: () => api("/api/view", { patch: { toggle: { group: c.group, key: c.key } } }),
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
  const series = currentSeries();
  const q = (s) => `"${s.replace(/"/g, '""')}"`;
  const T = result.axes.temperature_c, P = result.axes.pressure_bar;
  const rows = [[q(`${T.label} (${T.unit})`), q(`${P.label} (${P.unit})`), ...series.map((s) => q(s.csv))].join(",")];
  T.values.forEach((t, i) => rows.push([t, P.values[i], ...series.map((s) => (s.y[i] === null ? "" : s.y[i]))].join(",")));
  download(`${result.id}.csv`, rows.join("\n"), "text/csv");
});
$("dl-svg").addEventListener("click", () => download(`${result.id}.svg`, plot.svgText(), "image/svg+xml"));
