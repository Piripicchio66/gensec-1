/* ==========================================================
 * GenSec GUI — data store
 *
 * Backend-backed replacement of the former mock.  Exposes the
 * same window.GS_DATA shape the rest of the UI already uses
 * (SECTION, MATERIALS, DEMANDS, COMBINATIONS, ENVELOPES,
 *  VERIFICATION, makeMxMyContour, makeNM, makeMchi,
 *  makeSurfaceSlices), so no change is needed in plots.jsx /
 * panels.jsx / app.jsx.
 *
 * Two modes:
 *   - EMPTY (no YAML loaded yet): tables are empty, plots show
 *     a "Load a YAML to begin" placeholder, the Run button in
 *     the drawer becomes the primary CTA.
 *   - LOADED: GS_DATA is populated from the AnalysisResult
 *     returned by /api/analyze, and the curve generators pull
 *     from the numeric domain payload the backend shipped.
 *
 * The store is plain-window (no React context): plots.jsx reads
 * from window.GS_DATA synchronously.  App.jsx calls GS_STORE.set()
 * after a successful analyze(), which swaps GS_DATA and forces a
 * React re-render by bumping a generation counter.
 * ========================================================== */
(function () {
  "use strict";

  // ---- Empty defaults (used before first successful analyze) ----
  const EMPTY = {
    SECTION: { B: 0, H: 0, bulk_material: "—",
               n_fibers_x: 0, n_fibers_y: 0, rebars: [] },
    MATERIALS: [],
    DEMANDS: [],
    COMBINATIONS: [],
    ENVELOPES: [],
    VERIFICATION: [],
    // Curve generators return empty arrays in the empty state.
    makeMxMyContour: () => [],
    makeNM: () => [],
    makeMchi: () => [],
    makeSurfaceSlices: () => [],
    // Extra payload used by the plot tab captions.
    _loaded: false,
  };

  // ---- Conversion helpers ----

  // Safely read a field on a Pydantic-like nested object.
  const g = (o, k, dflt) => (o && o[k] !== undefined && o[k] !== null) ? o[k] : dflt;

  /** Build Mx-My contour generator from the backend domain payload. */
  function makeContourFactory(domain) {
    // Backend may ship either (a) a single "nm_mx_my" list keyed by N,
    // or (b) nothing -- in which case we fall back to linear scaling
    // of the resistance-only NM diagram (rough but non-empty).
    const slices = g(domain, "mx_my_slices", null);  // [{N_kN, points:[[Mx,My],...]}]
    if (slices && slices.length) {
      // Return nearest-N slice; perfect for a slider.
      const sorted = [...slices].sort((a, b) => a.N_kN - b.N_kN);
      return function (N_kN /*, n */) {
        let best = sorted[0], bd = Math.abs(N_kN - best.N_kN);
        for (const s of sorted) {
          const d = Math.abs(N_kN - s.N_kN);
          if (d < bd) { best = s; bd = d; }
        }
        return best.points && best.points.length
          ? best.points.concat([best.points[0]])
          : [];
      };
    }
    // Fallback: scale the uniaxial rectangle envelope.
    const nm = g(domain, "nm_points", null);
    if (nm && nm.length) {
      return function (N_kN) {
        // Find the N closest to N_kN and use its M as the Mx bound,
        // then draw an ellipse Mx/My with My capped at 0.6*Mx.
        let Mx = 0;
        let bd = Infinity;
        for (const [N, M] of nm) {
          const d = Math.abs(N - N_kN);
          if (d < bd) { bd = d; Mx = Math.abs(M); }
        }
        const My = Mx * 0.6;
        const pts = [];
        const n = 144;
        for (let i = 0; i < n; i++) {
          const t = (i / n) * Math.PI * 2;
          pts.push([Mx * Math.cos(t), My * Math.sin(t)]);
        }
        pts.push(pts[0]);
        return pts;
      };
    }
    return () => [];
  }

  /** Uniaxial N–M closed contour from the nm_points list (upper branch only). */
  function makeNMFactory(domain) {
    const nm = g(domain, "nm_points", null);
    if (!nm || !nm.length) return () => [];
    return function () {
      const upper = nm.map(([N, M]) => [N, M]);
      const lower = [...upper].reverse().map(([N, M]) => [N, -M]);
      return upper.concat(lower).concat([upper[0]]);
    };
  }

  /** M-χ curves per N level (backend ships as array of {N, points}). */
  function makeMchiFactory(domain) {
    const series = g(domain, "mchi", null);
    if (!series || !series.length) return () => [];
    return function () {
      return series.map(s => ({
        N: s.N_kN,
        pts: (s.points || []).map(p => [p[0], p[1]]),
      }));
    };
  }

  /** 3D surface as stacked contours.  Lazy: may be null until the user opens the tab. */
  function makeSurfaceFactory(domain) {
    const slices = g(domain, "surface_slices", null);
    if (!slices || !slices.length) return () => [];
    return function () {
      return slices.map(s => ({
        N: s.N_kN,
        contour: (s.points || []).concat((s.points && s.points.length) ? [s.points[0]] : []),
      }));
    };
  }

  // ---- Project AnalysisResult -> GS_DATA shape ----

  function fromAnalysisResult(ar) {
    const section = ar.section || {};
    const domain = ar.domain || {};
    const verif = ar.verification || {};

    const SECTION = {
      B: g(section, "B_mm", 0),
      H: g(section, "H_mm", 0),
      bulk_material: g(section, "bulk_material", "—"),
      n_fibers_x: g(section, "n_fibers_x", 0),
      n_fibers_y: g(section, "n_fibers_y", 0),
      rebars: (section.rebars || []).map(r => ({
        x: r.x_mm, y: r.y_mm,
        diameter: r.diameter_mm, material: r.material,
        As: r.area_mm2,
      })),
    };

    const MATERIALS = (ar.materials || []).map(m => {
      const props = m.properties || {};
      return {
        id: m.id,
        kind: m.kind,
        cls: g(props, "class", m.id),
        fcd: props.fcd_MPa ? props.fcd_MPa.toFixed(2) + " MPa" : undefined,
        fctm: props.fctm_MPa ? props.fctm_MPa.toFixed(2) + " MPa" : undefined,
        Ecm: props.Ecm_MPa ? props.Ecm_MPa.toLocaleString() + " MPa" : undefined,
        fyd: props.fyd_MPa ? props.fyd_MPa.toFixed(1) + " MPa" : undefined,
        Es:  props.Es_MPa  ? props.Es_MPa.toLocaleString() + " MPa" : undefined,
        eps_su: props.eps_su ? props.eps_su.toString() : undefined,
      };
    });

    const DEMANDS = (ar.demands || []).map(d => ({
      name: d.name, N: d.N_kN, Mx: d.Mx_kNm, My: d.My_kNm,
    }));

    const COMBINATIONS = (ar.combinations || []).map(c => ({
      name: c.name,
      staged: !!c.staged,
      components: c.components || [],
      stages: c.stages || [],
      resolved: {
        N: g(c.resolved, "N_kN", 0),
        Mx: g(c.resolved, "Mx_kNm", 0),
        My: g(c.resolved, "My_kNm", 0),
      },
    }));

    const ENVELOPES = (ar.envelopes || []).map(e => ({
      name: e.name,
      members: e.members || [],
      eta_max: e.eta_max,
    }));

    // Verification rows: backend already emits flat rows mixing
    // demand/combination/envelope.  Ensure display fields are set.
    const VERIFICATION = (verif.rows || []).map(r => {
      const N = (r.N_kN !== null && r.N_kN !== undefined) ? r.N_kN : "—";
      const Mx = (r.Mx_kNm !== null && r.Mx_kNm !== undefined) ? r.Mx_kNm : "—";
      const My = (r.My_kNm !== null && r.My_kNm !== undefined) ? r.My_kNm : "—";
      return {
        kind: r.kind, name: r.name,
        N, Mx, My,
        eta3D: r.eta_3D, eta2D: r.eta_2D,
        etaPath: r.eta_path, etaPath2D: r.eta_path_2D,
        status: r.status || "ok",
        staged: !!r.staged,
      };
    });

    return {
      SECTION, MATERIALS, DEMANDS, COMBINATIONS, ENVELOPES, VERIFICATION,
      makeMxMyContour: makeContourFactory(domain),
      makeNM:          makeNMFactory(domain),
      makeMchi:        makeMchiFactory(domain),
      makeSurfaceSlices: makeSurfaceFactory(domain),
      _loaded: true,
      _meta: ar.meta || {},
    };
  }

  // ---- Public store ----

  const listeners = new Set();

  const GS_STORE = {
    isLoaded() { return !!(window.GS_DATA && window.GS_DATA._loaded); },
    setFromAnalysis(ar) {
      window.GS_DATA = fromAnalysisResult(ar);
      listeners.forEach(fn => { try { fn(); } catch (_) {} });
    },
    clear() {
      window.GS_DATA = Object.assign({}, EMPTY);
      listeners.forEach(fn => { try { fn(); } catch (_) {} });
    },
    subscribe(fn) { listeners.add(fn); return () => listeners.delete(fn); },
  };

  window.GS_DATA = Object.assign({}, EMPTY);
  window.GS_STORE = GS_STORE;
})();
