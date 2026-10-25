"use strict";

/* A small line plot: linear or log y axis, gaps for missing points, hover readout, SVG export.
   Series share one x array unless they carry their own (a compared run).
   series = [{label, color, dash, y: [...], x?: [...], lo?: [...], hi?: [...]}]; lo and hi shade a range. */
class LinePlot {
  constructor(host) {
    this.host = host;
    this.data = null;
    this.tip = h("div", { class: "tip" });
    this.host.append(this.tip);
    this.svgHost = h("div", { style: "position:absolute;inset:0" });
    this.host.prepend(this.svgHost);
    this.geom = null;
    new ResizeObserver(() => this.render()).observe(host);
    host.addEventListener("mousemove", (e) => this.hover(e));
    host.addEventListener("mouseleave", () => this.leave());
  }

  set(x, series, opts) {
    this.data = { x, series, opts: opts || {} };
    this.render();
  }

  static niceTicks(lo, hi, n = 6) {
    const span = hi - lo || Math.abs(hi) || 1;
    const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
    const first = Math.ceil(lo / step - 1e-9) * step, out = [];
    for (let v = first; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    return out;
  }

  render() {
    const d = this.data;
    const w = this.host.clientWidth, hgt = this.host.clientHeight;
    if (!d || w < 50 || hgt < 50) return;
    const { x, series, opts } = d;
    const m = { l: 70, r: 18, t: 12, b: 46 };
    const pw = w - m.l - m.r, ph = hgt - m.t - m.b;
    const xs = [x, ...series.map((s) => s.x || [])].flat().filter((v) => v !== null && isFinite(v));
    const x0 = Math.min(...xs), x1 = Math.max(...xs);

    let ys = [];
    const keep = (v) => { if (v !== null && isFinite(v) && (!opts.log || v > 0)) ys.push(v); };
    series.forEach((s) => [s.y, s.lo || [], s.hi || []].forEach((arr) => arr.forEach(keep)));
    if (opts.refLines) ys = ys.concat(opts.log ? [] : opts.refLines);
    let y0 = ys.length ? Math.min(...ys) : 0, y1 = ys.length ? Math.max(...ys) : 1;
    let yticks;
    if (opts.log) {
      let e0 = Math.floor(Math.log10(y0)), e1 = Math.ceil(Math.log10(y1));
      if (e1 === e0) e1 += 1;
      if (e1 - e0 > 14) e0 = e1 - 14;
      y0 = Math.pow(10, e0); y1 = Math.pow(10, e1);
      yticks = [];
      const stride = Math.ceil((e1 - e0) / 8);
      for (let e = e0; e <= e1; e += stride) yticks.push(Math.pow(10, e));
    } else {
      if (y0 >= 0 && y0 < 0.35 * y1) y0 = 0;
      if (y1 === y0) { y1 += 1; y0 -= 1; }
      const pad = (y1 - y0) * 0.05;
      if (y0 !== 0) y0 -= pad;
      y1 += pad;
      yticks = LinePlot.niceTicks(y0, y1);
    }
    const sx = (v) => m.l + ((v - x0) / (x1 - x0 || 1)) * pw;
    const sy = opts.log
      ? (v) => m.t + ph - ((Math.log10(v) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0))) * ph
      : (v) => m.t + ph - ((v - y0) / (y1 - y0)) * ph;
    const xticks = LinePlot.niceTicks(x0, x1, 8);

    const f = (n) => n.toFixed(1);
    let g = "";
    yticks.forEach((t) => {
      const yy = sy(t);
      if (yy < m.t - 1 || yy > m.t + ph + 1) return;
      g += `<line class="gridline" x1="${m.l}" x2="${m.l + pw}" y1="${f(yy)}" y2="${f(yy)}"/>`;
      const label = opts.log ? `10${sup(Math.round(Math.log10(t)))}` : fmtVal(t);
      g += `<text x="${m.l - 8}" y="${f(yy + 4)}" text-anchor="end">${label}</text>`;
    });
    xticks.forEach((t) => {
      g += `<line class="axis" x1="${f(sx(t))}" x2="${f(sx(t))}" y1="${m.t + ph}" y2="${m.t + ph + 4}"/>`;
      g += `<text x="${f(sx(t))}" y="${m.t + ph + 18}" text-anchor="middle">${fmtVal(t)}</text>`;
    });
    (opts.refLines || []).forEach((v) => {
      if (!opts.log) g += `<line class="refline" x1="${m.l}" x2="${m.l + pw}" y1="${f(sy(v))}" y2="${f(sy(v))}"/>`;
    });
    // reliability strip along the bottom edge: the colour of each point's stretch of the path
    if (opts.strip) {
      const fills = { good: "var(--ok)", caution: "var(--warn)", poor: "var(--bad)" };
      const pts = x.map((v, i) => ({ v, i })).filter((p) => p.v !== null);
      pts.forEach((p, k) => {
        const lo = k === 0 ? p.v : (pts[k - 1].v + p.v) / 2, hi = k === pts.length - 1 ? p.v : (p.v + pts[k + 1].v) / 2;
        const xa = Math.min(sx(lo), sx(hi)), wd = Math.max(Math.abs(sx(hi) - sx(lo)), 2);
        g += `<rect x="${f(xa)}" y="${m.t + ph - 5}" width="${f(wd)}" height="5" ` +
          `style="fill:${fills[opts.strip[p.i]] || fills.good};stroke:none" opacity="0.85"/>`;
      });
    }
    // uncertainty ranges under every line: one filled band per unbroken stretch of the path
    series.forEach((s) => {
      if (!s.lo || !s.hi) return;
      const sxs = s.x || x;
      const ok = (i) => s.lo[i] !== null && s.hi[i] !== null && sxs[i] !== null && sxs[i] !== undefined &&
        isFinite(s.lo[i]) && isFinite(s.hi[i]) && (!opts.log || s.hi[i] > 0);
      const clamp = (v) => Math.min(Math.max(v, y0), y1);           // a zero low edge sits on the axis
      let run = [];
      const flush = () => {
        if (run.length > 1) {
          const top = run.map((i) => `${f(sx(sxs[i]))},${f(sy(clamp(s.hi[i])))}`);
          const bot = run.slice().reverse().map((i) => `${f(sx(sxs[i]))},${f(sy(clamp(s.lo[i])))}`);
          g += `<path class="band" d="M${top.join("L")}L${bot.join("L")}Z" style="fill:${s.color};stroke:none" ` +
            `fill-opacity="${s.dash ? 0.12 : 0.2}"/>`;
        }
        run = [];
      };
      s.y.forEach((_, i) => (ok(i) ? run.push(i) : flush()));
      flush();
    });
    // thick solid lines first, then thin darker dashed ones on top, so a dashed line that follows
    // a solid one exactly is still visible
    [...series].sort((a, b) => (a.dash ? 1 : 0) - (b.dash ? 1 : 0)).forEach((s) => {
      let path = "", pen = false;
      const sxs = s.x || x;
      s.y.forEach((v, i) => {
        if (v === null || !isFinite(v) || (opts.log && v <= 0) || sxs[i] === null || sxs[i] === undefined) { pen = false; return; }
        path += `${pen ? "L" : "M"}${f(sx(sxs[i]))},${f(sy(v))}`;
        pen = true;
      });
      const stroke = s.dash ? `color-mix(in srgb, ${s.color} 55%, var(--ink))` : s.color;
      g += `<path d="${path}" fill="none" style="stroke:${stroke}" stroke-width="${s.dash ? 1.7 : 3.2}" ` +
        `stroke-linejoin="round" stroke-linecap="round"${s.dash ? ` stroke-dasharray="${s.dash}"` : ""}/>`;
    });
    this.svgHost.innerHTML =
      `<svg viewBox="0 0 ${w} ${hgt}" width="${w}" height="${hgt}" xmlns="http://www.w3.org/2000/svg" role="img" ` +
      `aria-label="${(opts.title || "Plot").replace(/"/g, "")}">` +
      `<g class="axis"><line x1="${m.l}" x2="${m.l + pw}" y1="${m.t + ph}" y2="${m.t + ph}"/>` +
      `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ph}"/>${g}</g>` +
      `<text class="axislabel" x="${m.l + pw / 2}" y="${hgt - 6}" text-anchor="middle">${opts.xLabel || ""}</text>` +
      `<text class="axislabel" transform="translate(14,${m.t + ph / 2}) rotate(-90)" text-anchor="middle">${opts.yLabel || ""}</text>` +
      `<line class="cross" id="cross" y1="${m.t}" y2="${m.t + ph}" style="display:none"/></svg>`;
    this.geom = { m, pw, ph, x0, x1, sx, sy };
  }

  hover(e) {
    if (!this.data || !this.geom) return;
    const { x, series, opts } = this.data, { m, pw, sx, sy } = this.geom;
    const rect = this.host.getBoundingClientRect();
    const px = e.clientX - rect.left;
    if (px < m.l || px > m.l + pw) return this.leave();
    let best = -1, bd = Infinity;
    x.forEach((v, i) => { if (v !== null && Math.abs(sx(v) - px) < bd) { bd = Math.abs(sx(v) - px); best = i; } });
    if (best < 0) return;
    const cross = this.svgHost.querySelector("#cross");
    if (cross) { cross.setAttribute("x1", sx(x[best])); cross.setAttribute("x2", sx(x[best])); cross.style.display = ""; }
    let html = `<b>${opts.xName || "x"} ${fmtVal(x[best])}${opts.xUnit || ""}</b>`;
    if (opts.sub) html += `<div class="muted">${opts.sub(best)}</div>`;
    if (opts.strip && opts.strip[best] && opts.strip[best] !== "good") {
      html += `<div class="muted">thermodynamics: ${opts.strip[best] === "poor" ? "unreliable" : "approximate"}</div>`;
    }
    series.forEach((s) => {
      let k = best;
      if (s.x) {          // its own grid: the point nearest the cursor's x
        let d = Infinity;
        s.x.forEach((v, i) => { if (v !== null && Math.abs(v - x[best]) < d) { d = Math.abs(v - x[best]); k = i; } });
      }
      const range = s.lo && s.lo[k] !== null && s.hi[k] !== null ? ` <span class="muted">(${fmtVal(s.lo[k])} to ${fmtVal(s.hi[k])})</span>` : "";
      html += `<div><span style="color:${s.color}">■</span> ${s.label}: ${fmtVal(s.y[k])}${range}</div>`;
    });
    this.tip.innerHTML = html;
    this.tip.style.display = "block";
    const tw = this.tip.offsetWidth;
    this.tip.style.left = `${Math.min(Math.max(sx(x[best]) + 12, 4), rect.width - tw - 4)}px`;
    this.tip.style.top = `${m.t + 6}px`;
  }

  leave() {
    this.tip.style.display = "none";
    const cross = this.svgHost.querySelector("#cross");
    if (cross) cross.style.display = "none";
  }

  /* A standalone copy of the drawing, with the page's colours written into it. */
  svgText() {
    const svg = this.svgHost.querySelector("svg");
    if (!svg) return "";
    const clone = svg.cloneNode(true);
    const src = svg.querySelectorAll("*"), dst = clone.querySelectorAll("*");
    src.forEach((el, i) => {
      const cs = getComputedStyle(el);
      let style = "";
      for (const p of ["stroke", "fill", "stroke-width", "stroke-dasharray", "font-size", "font-weight",
        "font-style", "opacity", "fill-opacity", "paint-order", "stroke-linejoin"]) {
        style += `${p}:${cs.getPropertyValue(p)};`;
      }
      style += `font-family:${cs.fontFamily};`;
      dst[i].setAttribute("style", style);
    });
    const bg = getComputedStyle(document.body).backgroundColor;
    clone.insertAdjacentHTML("afterbegin", `<rect width="100%" height="100%" fill="${bg}"/>`);
    return new XMLSerializer().serializeToString(clone);
  }
}
