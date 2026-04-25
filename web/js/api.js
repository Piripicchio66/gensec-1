/* ==========================================================
 * GenSec GUI — HTTP client
 *
 * Thin wrapper over fetch() that targets the FastAPI backend
 * (gensec.server).  Every call posts a JSON body, parses the
 * JSON response, and throws a friendly Error on non-2xx.
 *
 * All endpoints share the same convention: the YAML source is
 * passed as a plain string in `yaml_text` -- never as a path.
 * That keeps the server stateless and lets us send unsaved
 * buffers straight from the inline editor.
 * ========================================================== */
(function () {
  "use strict";

  // Relative to the page origin -- works in dev (same uvicorn) and in
  // a frozen desktop build identically.  Override via ?api=... for
  // manual testing (e.g. hitting a remote instance from file://).
  const urlParams = new URLSearchParams(location.search);
  const BASE = (urlParams.get("api") || "").replace(/\/$/, "");

  async function postJSON(path, body) {
    const r = await fetch(BASE + path, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body || {}),
    });
    const text = await r.text();
    let data = null;
    try { data = text ? JSON.parse(text) : null; } catch (_) { /* non-JSON */ }
    if (!r.ok) {
      const detail = (data && data.detail) || text || r.statusText;
      const err = new Error(detail);
      err.status = r.status;
      throw err;
    }
    return data;
  }

  async function getJSON(path) {
    const r = await fetch(BASE + path);
    if (!r.ok) throw new Error(r.statusText);
    return r.json();
  }

  const GensecAPI = {
    /** Health / cache stats.  Returns {status, gensec_version, cache}. */
    health: () => getJSON("/api/health"),

    /** Run the full analysis.  Returns an AnalysisResult. */
    analyze: (yaml_text) => postJSON("/api/analyze", { yaml_text }),

    /** Mx–My contour at fixed N.  Cheap after first analyze(). */
    contour: (yaml_text, N_kN, n_angles = 144) =>
      postJSON("/api/contour", { yaml_text, N_kN, n_angles }),

    /** Verify one ad-hoc demand point. */
    verifyPoint: (yaml_text, N_kN, Mx_kNm, My_kNm = 0) =>
      postJSON("/api/verify-point", { yaml_text, N_kN, Mx_kNm, My_kNm }),

    /** Render a matplotlib plot as base64 PNG. */
    plot: (kind, yaml_text, options = {}, width_px = 1200, height_px = 800, dpi = 150) =>
      postJSON("/api/plot/" + encodeURIComponent(kind), {
        yaml_text, width_px, height_px, dpi, options,
      }),
  };

  window.GensecAPI = GensecAPI;
})();
