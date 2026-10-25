"use strict";

/* A predominance or activity diagram: fields from a grid of winners, boundaries along the cell edges,
   a hatched band where the winner's lead is within the logK uncertainty, dashed water and saturation
   lines, and the run's path on top.  data is the server's /api/diagram answer. */
class DiagramPlot {
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
    host.addEventListener("mouseleave", () => { this.tip.style.display = "none"; });
  }

  set(data) { this.data = data; this.render(); }

  render() {
    const d = this.data;
    const w = this.host.clientWidth, hgt = this.host.clientHeight;
    if (!d || w < 50 || hgt < 50) { this.svgHost.innerHTML = ""; return; }
    const m = { l: 64, r: 18, t: 12, b: 46 };
    const pw = w - m.l - m.r, ph = hgt - m.t - m.b;
    const { nx, ny, x, y } = d;
    const sx = (v) => m.l + ((v - x.lo) / (x.hi - x.lo)) * pw;
    const sy = (v) => m.t + ph - ((v - y.lo) / (y.hi - y.lo)) * ph;
    const cw = pw / (nx - 1), chh = ph / (ny - 1);
    // cell (i, j) is centred on grid point j along x and i along y
    const cx = (j) => Math.max(m.l, m.l + (j - 0.5) * cw), cx2 = (j) => Math.min(m.l + pw, m.l + (j + 0.5) * cw);
    const cy = (i) => Math.min(m.t + ph, m.t + ph - (i - 0.5) * chh), cy2 = (i) => Math.max(m.t, m.t + ph - (i + 0.5) * chh);
    const f = (n) => n.toFixed(1);
    const order = d.fields.map((fl) => fl.index);
    const colour = {};
    d.fields.forEach((fl, k) => {
      const share = fl.kind === "mineral" ? 30 : 14;
      colour[fl.index] = `color-mix(in srgb, var(--s${(k % 8) + 1}) ${share}%, var(--panel))`;
    });
    let g = "";
    // fields, as runs of equal cells along each row
    for (let i = 0; i < ny; i++) {
      const row = d.grid[i];
      let start = 0;
      for (let j = 1; j <= nx; j++) {
        if (j === nx || row[j] !== row[start]) {
          g += `<rect x="${f(cx(start))}" y="${f(cy2(i))}" width="${f(cx2(j - 1) - cx(start) + 0.4)}" ` +
            `height="${f(cy(i) - cy2(i) + 0.4)}" style="fill:${colour[row[start]]};stroke:none"/>`;
          start = j;
        }
      }
    }
    // the band where the boundary is uncertain at this (T, P)
    for (let i = 0; i < ny; i++) {
      const row = d.uncertain[i];
      let j = 0;
      while (j < nx) {
        if (row[j] !== "1") { j++; continue; }
        let k = j;
        while (k + 1 < nx && row[k + 1] === "1") k++;
        g += `<rect x="${f(cx(j))}" y="${f(cy2(i))}" width="${f(cx2(k) - cx(j) + 0.4)}" ` +
          `height="${f(cy(i) - cy2(i) + 0.4)}" fill="url(#band)"/>`;
        j = k + 1;
      }
    }
    // extrapolated fields: dotted
    const flagged = new Set(d.fields.filter((fl) => fl.extrapolated).map((fl) => fl.index));
    if (flagged.size) {
      for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
        if (flagged.has(d.grid[i][j])) g += `<rect x="${f(cx(j))}" y="${f(cy2(i))}" width="${f(cx2(j) - cx(j))}" height="${f(cy(i) - cy2(i))}" fill="url(#dots)"/>`;
      }
    }
    // boundaries along the cell edges
    let edges = "";
    for (let i = 0; i < ny; i++) for (let j = 0; j < nx; j++) {
      const v = d.grid[i][j];
      if (j + 1 < nx && d.grid[i][j + 1] !== v) edges += `M${f(cx2(j))},${f(cy(i))}V${f(cy2(i))}`;
      if (i + 1 < ny && d.grid[i + 1][j] !== v) edges += `M${f(cx(j))},${f(cy2(i))}H${f(cx2(j))}`;
    }
    g += `<path d="${edges}" fill="none" style="stroke:var(--ink)" stroke-width="1.1" opacity="0.75"/>`;
    // water stability and saturation limits
    const clip = `clip-path="url(#plotclip)"`;
    Object.entries(d.water || {}).forEach(([key, pts]) => {
      g += `<line ${clip} x1="${f(sx(pts[0][0]))}" y1="${f(sy(pts[0][1]))}" x2="${f(sx(pts[1][0]))}" y2="${f(sy(pts[1][1]))}" ` +
        `style="stroke:var(--accent)" stroke-width="1.6" stroke-dasharray="6 4"/>`;
    });
    (d.limits || []).forEach((lim) => {
      if (lim.axis === "x") {
        const xx = sx(lim.value);
        g += `<rect x="${f(xx)}" y="${m.t}" width="${f(m.l + pw - xx)}" height="${ph}" style="fill:var(--ink)" opacity="0.07"/>`;
        g += `<line x1="${f(xx)}" x2="${f(xx)}" y1="${m.t}" y2="${m.t + ph}" style="stroke:var(--accent)" stroke-width="1.6" stroke-dasharray="6 4"/>`;
        g += `<text class="small" x="${f(xx + 4)}" y="${m.t + 14}" style="fill:var(--accent)">${lim.label}</text>`;
      } else {
        const yy = sy(lim.value);
        g += `<rect x="${m.l}" y="${m.t}" width="${pw}" height="${f(yy - m.t)}" style="fill:var(--ink)" opacity="0.07"/>`;
        g += `<line x1="${m.l}" x2="${m.l + pw}" y1="${f(yy)}" y2="${f(yy)}" style="stroke:var(--accent)" stroke-width="1.6" stroke-dasharray="6 4"/>`;
        g += `<text class="small" x="${m.l + 6}" y="${f(yy - 5)}" style="fill:var(--accent)">${lim.label}</text>`;
      }
    });
    // labels
    d.fields.forEach((fl) => {
      if (fl.small) return;
      const [lx, ly] = fl.at;
      const style = fl.kind === "mineral" ? "font-weight:600" : "font-style:italic";
      g += `<text x="${f(sx(lx))}" y="${f(sy(ly) + 4)}" text-anchor="middle" class="flabel" style="${style}">` +
        `${esc(fl.label)}${fl.extrapolated ? "*" : ""}</text>`;
    });
    // the run's path
    Object.entries(d.path || {}).forEach(([mode, pts]) => {
      const line = pts.map((p, k) => `${k ? "L" : "M"}${f(sx(p[0]))},${f(sy(p[1]))}`).join("");
      g += `<path ${clip} d="${line}" fill="none" style="stroke:var(--ink)" stroke-width="2.2" ` +
        `${mode === "flow_through" ? 'stroke-dasharray="7 5"' : ""}/>`;
      pts.forEach((p, k) => {
        const near = Math.abs(p[2] - d.t_c) <= 10;
        g += `<circle ${clip} cx="${f(sx(p[0]))}" cy="${f(sy(p[1]))}" r="${k === pts.length - 1 ? 5 : 3.6}" ` +
          `style="fill:${near ? "var(--ink)" : "var(--panel)"};stroke:var(--ink)" stroke-width="1.5"/>`;
      });
    });
    // axes
    let ax = "";
    LinePlot.niceTicks(y.lo, y.hi).forEach((t) => {
      if (t < y.lo - 1e-9 || t > y.hi + 1e-9) return;
      ax += `<line x1="${m.l - 4}" x2="${m.l}" y1="${f(sy(t))}" y2="${f(sy(t))}"/>`;
      ax += `<text x="${m.l - 8}" y="${f(sy(t) + 4)}" text-anchor="end">${fmtVal(t)}</text>`;
    });
    LinePlot.niceTicks(x.lo, x.hi, 8).forEach((t) => {
      if (t < x.lo - 1e-9 || t > x.hi + 1e-9) return;
      ax += `<line x1="${f(sx(t))}" x2="${f(sx(t))}" y1="${m.t + ph}" y2="${m.t + ph + 4}"/>`;
      ax += `<text x="${f(sx(t))}" y="${m.t + ph + 18}" text-anchor="middle">${fmtVal(t)}</text>`;
    });
    this.svgHost.innerHTML =
      `<svg viewBox="0 0 ${w} ${hgt}" width="${w}" height="${hgt}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="${esc(d.title)}">` +
      `<defs><pattern id="band" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">` +
      `<line x1="0" y1="0" x2="0" y2="6" style="stroke:var(--warn)" stroke-width="1.6" opacity="0.4"/></pattern>` +
      `<pattern id="dots" width="5" height="5" patternUnits="userSpaceOnUse"><circle cx="2.5" cy="2.5" r="0.9" style="fill:var(--muted)"/></pattern>` +
      `<clipPath id="plotclip"><rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}"/></clipPath></defs>` +
      `<g>${g}</g>` +
      `<g class="axis"><rect x="${m.l}" y="${m.t}" width="${pw}" height="${ph}" fill="none"/>${ax}</g>` +
      `<text class="axislabel" x="${m.l + pw / 2}" y="${hgt - 6}" text-anchor="middle">${esc(x.label)}</text>` +
      `<text class="axislabel" transform="translate(14,${m.t + ph / 2}) rotate(-90)" text-anchor="middle">${esc(y.label)}</text></svg>`;
    this.geom = { m, pw, ph };
    this.order = order;
  }

  hover(e) {
    const d = this.data, g = this.geom;
    if (!d || !g) return;
    const r = this.host.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const fx = (px - g.m.l) / g.pw, fy = 1 - (py - g.m.t) / g.ph;
    if (fx < 0 || fx > 1 || fy < 0 || fy > 1) { this.tip.style.display = "none"; return; }
    const j = Math.round(fx * (d.nx - 1)), i = Math.round(fy * (d.ny - 1));
    const field = d.fields.find((fl) => fl.index === d.grid[i][j]);
    const xv = d.x.lo + fx * (d.x.hi - d.x.lo), yv = d.y.lo + fy * (d.y.hi - d.y.lo);
    const close = d.uncertain[i][j] === "1" ? "<div class='muted'>within the uncertainty of the boundary</div>" : "";
    this.tip.innerHTML = `<b>${esc(field ? field.label : "")}</b>${field && field.extrapolated ? " (data extrapolated)" : ""}` +
      `<div>${esc(d.x.label)} ${xv.toFixed(2)}, ${esc(d.y.label)} ${yv.toFixed(2)}</div>${close}`;
    this.tip.style.display = "block";
    this.tip.style.left = Math.min(px + 14, r.width - 220) + "px";
    this.tip.style.top = Math.max(py - 10, 0) + "px";
  }

}
DiagramPlot.prototype.svgText = LinePlot.prototype.svgText;

function esc(s) {
  return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}
