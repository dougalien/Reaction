"use strict";

async function api(path, body) {
  const r = await fetch(path, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body || {}),
  });
  return r.json();
}

function connect(onState, onStatus) {
  let es = null, retry = null;
  function open() {
    es = new EventSource("/events");
    es.onmessage = (e) => { onStatus && onStatus(true); onState(JSON.parse(e.data)); };
    es.onerror = () => {
      onStatus && onStatus(false);
      es.close();
      clearTimeout(retry);
      retry = setTimeout(open, 2000);
    };
  }
  open();
}

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (k === "class") e.className = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (v !== false && v !== null && v !== undefined) e.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat()) {
    if (kid === null || kid === undefined || kid === false) continue;
    e.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return e;
}

function setIfIdle(input, value) {
  if (document.activeElement !== input) input.value = value;
}

const SUP = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵",
  "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻" };
function sup(n) { return String(n).split("").map((c) => SUP[c] || c).join(""); }

function fmtVal(v) {
  if (v === null || v === undefined || !isFinite(v)) return "–";
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a < 1e-3 || a >= 1e5) {
    const [m, e] = v.toExponential(2).split("e");
    return `${m}×10${sup(parseInt(e, 10))}`;
  }
  return Number(v.toPrecision(3)).toString();
}
