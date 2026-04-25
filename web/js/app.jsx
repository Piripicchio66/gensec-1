/* ==========================================================
 * GenSec GUI — main app (backend-connected)
 *
 * Differences from the mock version:
 *   - window.GS_DATA starts EMPTY; populated only after the user
 *     loads a YAML (via file picker, drag-drop, load-example, or
 *     localStorage recall).
 *   - The Run button POSTs the current YAML buffer to /api/analyze
 *     and calls GS_STORE.setFromAnalysis() on success.
 *   - Errors go to a toast banner; loading state goes to the topbar
 *     dot + the empty-state card.
 *   - localStorage remembers the last YAML text + path so a reload
 *     reopens the session.
 * ========================================================== */

const { useState, useEffect, useMemo, useRef, useCallback } = React;

function useSize(ref) {
  const [size, setSize] = useState({ w: 320, h: 240 });
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((entries) => {
      for (const e of entries) {
        const cr = e.contentRect;
        setSize({ w: Math.max(100, cr.width), h: Math.max(100, cr.height) });
      }
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, [ref]);
  return size;
}

const DEFAULT_FLAGS = {
  eta_3D: true, eta_2D: true, eta_path: true, eta_path_2D: true,
  delta_N_tol: 0.03,
  generate_mx_my: true, generate_3d_surface: true,
  generate_moment_curvature: true, generate_polar_ductility: true,
  generate_3d_moment_curvature: true,
  n_angles_mx_my: 144,
};

const TWEAK_DEFAULTS = /*EDITMODE-BEGIN*/{
  "accent_h": 225,
  "dark": false,
  "serif": true,
  "compact": false,
  "grid": true,
  "rail_w": 260
}/*EDITMODE-END*/;

const LS_YAML_TEXT = "gensec.yamlText";
const LS_YAML_NAME = "gensec.yamlName";

// ----- Empty-state / drop-zone panel -----
function EmptyState({ onLoadFile, onLoadExample, loading, error }) {
  const fileInputRef = useRef(null);
  const [dragOver, setDragOver] = useState(false);

  return (
    <div style={{
      position: "absolute", inset: 0,
      display: "grid", placeItems: "center",
      background: "var(--paper)",
    }}>
      <div
        onDragOver={(e) => { e.preventDefault(); setDragOver(true); }}
        onDragLeave={() => setDragOver(false)}
        onDrop={(e) => {
          e.preventDefault();
          setDragOver(false);
          const f = e.dataTransfer.files && e.dataTransfer.files[0];
          if (f) onLoadFile(f);
        }}
        onClick={() => fileInputRef.current && fileInputRef.current.click()}
        style={{
          width: "min(520px, 80%)", padding: "48px 40px",
          border: `2px dashed ${dragOver ? "var(--accent)" : "var(--rule-2)"}`,
          borderRadius: "12px",
          background: dragOver ? "var(--accent-soft)" : "var(--paper-2)",
          textAlign: "center",
          cursor: "pointer",
          transition: "all 0.15s",
        }}
      >
        <input ref={fileInputRef} type="file" accept=".yaml,.yml"
               style={{ display: "none" }}
               onChange={(e) => e.target.files[0] && onLoadFile(e.target.files[0])} />
        <div style={{
          fontFamily: "var(--ff-serif)", fontSize: 22, fontWeight: 600,
          marginBottom: 10,
        }}>
          Load a GenSec YAML to begin
        </div>
        <div style={{ color: "var(--ink-3)", fontSize: 13, marginBottom: 20 }}>
          Drop a <code style={{ fontFamily: "var(--ff-mono)" }}>.yaml</code> file here, or click to browse
        </div>
        <button className="btn primary" onClick={(e) => { e.stopPropagation(); fileInputRef.current.click(); }}>
          Open YAML…
        </button>
        <span style={{ margin: "0 12px", color: "var(--ink-3)" }}>or</span>
        <button className="btn" onClick={(e) => { e.stopPropagation(); onLoadExample(); }}>
          Load example
        </button>
        {loading && (
          <div style={{ marginTop: 24, color: "var(--ink-3)", fontFamily: "var(--ff-mono)", fontSize: 12 }}>
            <span className="dot" style={{
              width: 8, height: 8, borderRadius: "50%",
              background: "var(--warn)", display: "inline-block", marginRight: 8,
            }} />
            Solving…
          </div>
        )}
        {error && (
          <div style={{
            marginTop: 24, padding: "10px 14px",
            background: "var(--fail-soft)", color: "var(--fail)",
            borderRadius: 6, fontFamily: "var(--ff-mono)", fontSize: 11.5,
            textAlign: "left", whiteSpace: "pre-wrap",
          }}>
            {error}
          </div>
        )}
      </div>
    </div>
  );
}

// ----- YAML editor modal -----
function YamlEditor({ open, initialText, initialName, onClose, onSubmit }) {
  const [text, setText] = useState(initialText || "");
  const [name, setName] = useState(initialName || "untitled.yaml");
  useEffect(() => {
    if (open) { setText(initialText || ""); setName(initialName || "untitled.yaml"); }
  }, [open, initialText, initialName]);

  if (!open) return null;
  return (
    <div style={{
      position: "fixed", inset: 0, background: "oklch(0% 0 0 / 0.35)",
      display: "grid", placeItems: "center", zIndex: 100,
    }} onClick={onClose}>
      <div onClick={(e) => e.stopPropagation()} style={{
        width: "min(900px, 92vw)", height: "min(80vh, 720px)",
        background: "var(--paper)", border: "1px solid var(--rule-2)",
        borderRadius: 10, boxShadow: "var(--sh-2)",
        display: "flex", flexDirection: "column", overflow: "hidden",
      }}>
        <div style={{
          display: "flex", alignItems: "center", gap: 10,
          padding: "10px 14px", borderBottom: "1px solid var(--rule)",
          background: "var(--paper-2)",
        }}>
          <h2 style={{ flex: 1 }}>Edit YAML</h2>
          <input value={name} onChange={(e) => setName(e.target.value)}
                 style={{
                   width: 240, height: 26, padding: "0 8px",
                   border: "1px solid var(--rule-2)", borderRadius: 4,
                   background: "var(--paper)",
                   fontFamily: "var(--ff-mono)", fontSize: 12,
                 }} />
          <button className="btn ghost sm" onClick={onClose}>Cancel</button>
          <button className="btn primary" onClick={() => onSubmit(text, name)}>
            Analyze
          </button>
        </div>
        <textarea
          value={text}
          onChange={(e) => setText(e.target.value)}
          spellCheck={false}
          style={{
            flex: 1, width: "100%", resize: "none",
            padding: 16, border: "none", outline: "none",
            background: "var(--paper)", color: "var(--ink)",
            fontFamily: "var(--ff-mono)", fontSize: 12.5, lineHeight: 1.55,
            tabSize: 2,
          }}
        />
      </div>
    </div>
  );
}

// ----- Error banner -----
function ErrorBanner({ message, onDismiss }) {
  if (!message) return null;
  return (
    <div style={{
      position: "fixed", top: 52, right: 16, zIndex: 90,
      maxWidth: 420, padding: "10px 14px",
      background: "var(--fail-soft)", color: "var(--fail)",
      border: "1px solid var(--fail)", borderRadius: 6,
      boxShadow: "var(--sh-2)",
      fontFamily: "var(--ff-mono)", fontSize: 11.5,
      display: "flex", alignItems: "flex-start", gap: 8,
    }}>
      <div style={{ flex: 1, whiteSpace: "pre-wrap" }}>{message}</div>
      <button className="btn ghost sm" onClick={onDismiss}
              style={{ color: "var(--fail)" }}>×</button>
    </div>
  );
}

function App() {
  // ---- Tabs / selection ----
  const [tab, setTab] = useState("section");
  const [resTab, setResTab] = useState("verify");
  const [selected, setSelected] = useState(null);
  const [flags, setFlags] = useState(DEFAULT_FLAGS);
  const [nSlice, setNSlice] = useState(0);

  // ---- Tweaks ----
  const [tweaks, setTweaks] = useState(TWEAK_DEFAULTS);
  const [tweaksOpen, setTweaksOpen] = useState(false);

  // ---- Backend state ----
  const [yamlText, setYamlText] = useState("");
  const [yamlName, setYamlName] = useState("(no file)");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [editorOpen, setEditorOpen] = useState(false);
  // Version counter bumped whenever GS_DATA changes, so React re-renders.
  const [dataVersion, setDataVersion] = useState(0);

  // Subscribe to store mutations.
  useEffect(() => window.GS_STORE.subscribe(() => setDataVersion(v => v + 1)), []);

  // Apply tweaks to CSS.
  useEffect(() => {
    document.documentElement.setAttribute("data-theme", tweaks.dark ? "dark" : "light");
    document.documentElement.style.setProperty("--accent-h", tweaks.accent_h);
    document.documentElement.style.setProperty("--rail-w", tweaks.rail_w + "px");
    document.documentElement.style.setProperty("--row-h", tweaks.compact ? "24px" : "28px");
    document.documentElement.style.setProperty("--pad", tweaks.compact ? "10px" : "14px");
    if (!tweaks.serif) document.documentElement.style.setProperty("--ff-serif", 'var(--ff-sans)');
    else document.documentElement.style.removeProperty("--ff-serif");
  }, [tweaks]);

  // Edit-mode protocol: register listener first, then announce.
  useEffect(() => {
    const handler = (ev) => {
      const d = ev.data;
      if (!d || typeof d !== "object") return;
      if (d.type === "__activate_edit_mode") setTweaksOpen(true);
      if (d.type === "__deactivate_edit_mode") setTweaksOpen(false);
    };
    window.addEventListener("message", handler);
    window.parent.postMessage({ type: "__edit_mode_available" }, "*");
    return () => window.removeEventListener("message", handler);
  }, []);

  const setTweak = (k, v) => {
    setTweaks((t) => {
      const next = { ...t, [k]: v };
      try {
        window.parent.postMessage({ type: "__edit_mode_set_keys", edits: { [k]: v } }, "*");
      } catch {}
      return next;
    });
  };
  const setFlag = (k, v) => setFlags((f) => ({ ...f, [k]: v }));

  // ---- Backend calls ----

  const runAnalyze = useCallback(async (text, name) => {
    if (!text || !text.trim()) {
      setError("YAML is empty. Load a file or paste content into the editor.");
      return;
    }
    setLoading(true);
    setError(null);
    try {
      const res = await window.GensecAPI.analyze(text);
      window.GS_STORE.setFromAnalysis(res);
      setYamlText(text);
      setYamlName(name || "(buffer)");
      try {
        localStorage.setItem(LS_YAML_TEXT, text);
        localStorage.setItem(LS_YAML_NAME, name || "(buffer)");
      } catch (_) {}
      // First demand becomes the selection; N-slice centers on it.
      const D = window.GS_DATA;
      if (D.DEMANDS.length) {
        setSelected("dem:" + D.DEMANDS[0].name);
        setNSlice(Math.round(D.DEMANDS[0].N || 0));
      } else {
        setSelected(null);
      }
    } catch (e) {
      setError(e.message || String(e));
    } finally {
      setLoading(false);
    }
  }, []);

  const handleLoadFile = useCallback(async (file) => {
    const text = await file.text();
    runAnalyze(text, file.name);
  }, [runAnalyze]);

  const handleLoadExample = useCallback(async () => {
    setLoading(true);
    setError(null);
    try {
      const r = await fetch("examples/biaxial_column.yaml");
      if (!r.ok) throw new Error("No example bundled (expected web/examples/biaxial_column.yaml).");
      const text = await r.text();
      await runAnalyze(text, "biaxial_column.yaml");
    } catch (e) {
      setError(e.message || String(e));
      setLoading(false);
    }
  }, [runAnalyze]);

  // Auto-load last YAML on first mount.
  useEffect(() => {
    const saved = (() => {
      try { return localStorage.getItem(LS_YAML_TEXT); } catch (_) { return null; }
    })();
    const savedName = (() => {
      try { return localStorage.getItem(LS_YAML_NAME) || "(buffer)"; } catch (_) { return "(buffer)"; }
    })();
    if (saved && saved.trim()) {
      setYamlText(saved);
      setYamlName(savedName);
      runAnalyze(saved, savedName);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // ---- Derived data ----
  const D = window.GS_DATA;
  // eslint-disable-next-line no-unused-vars
  const _ = dataVersion;  // force re-read when store bumps

  const isEmpty = !window.GS_STORE.isLoaded();

  const activeDemand = useMemo(() => {
    if (!selected) return null;
    if (selected.startsWith("dem:")) return D.DEMANDS.find((d) => d.name === selected.slice(4));
    if (selected.startsWith("cmb:")) {
      const c = D.COMBINATIONS.find((c) => c.name === selected.slice(4));
      if (c) return { name: c.name, ...c.resolved };
    }
    return null;
  }, [selected, dataVersion]);

  const activeName = activeDemand?.name;

  const demandPoints = useMemo(() => {
    const base = D.DEMANDS.map((d) => ({ name: d.name, N: d.N, Mx: d.Mx, My: d.My }));
    const combos = D.COMBINATIONS.map((c) => ({ name: c.name, ...c.resolved }));
    return base.concat(combos);
  }, [dataVersion]);

  const stageRef = useRef(null);
  const stageSize = useSize(stageRef);

  const tabs = [
    { id: "surface", label: "3D surface", cap: "N–Mx–My resistance surface" },
    { id: "mxmy",    label: "Mx–My",      cap: "Mx–My interaction contour at fixed N" },
    { id: "nm",      label: "N–M",        cap: "Uniaxial N–M interaction" },
    { id: "mchi",    label: "M–χ",        cap: "Moment-curvature" },
    { id: "polar",   label: "Polar",      cap: "Polar ductility rose" },
    { id: "section", label: "Section",    cap: "Fiber mesh & rebars" },
  ];

  const resTabs = [
    { id: "verify",       label: "Verification", count: D.VERIFICATION.length },
    { id: "demands",      label: "Demands",      count: D.DEMANDS.length },
    { id: "combinations", label: "Combinations", count: D.COMBINATIONS.length },
    { id: "envelopes",    label: "Envelopes",    count: D.ENVELOPES.length },
    { id: "materials",    label: "Materials",    count: D.MATERIALS.length },
  ];

  const resRows = useMemo(() => {
    if (resTab === "verify")       return D.VERIFICATION;
    if (resTab === "demands")      return D.VERIFICATION.filter(r => r.kind === "demand");
    if (resTab === "combinations") return D.VERIFICATION.filter(r => r.kind === "combo");
    if (resTab === "envelopes")    return D.VERIFICATION.filter(r => r.kind === "envelope");
    return [];
  }, [resTab, dataVersion]);

  const handleRun = () => runAnalyze(yamlText, yamlName);

  return (
    <div className="app" data-screen-label="GenSec GUI">
      {/* TOPBAR */}
      <div className="topbar">
        <div className="brand">
          <span className="logo">G</span>
          <h1>GenSec</h1>
          <span className="ver">v0.3</span>
        </div>
        <div className="topbar-sep"/>
        <div className="topbar-file">
          <svg className="ic" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.4">
            <path d="M3 2h7l3 3v9H3z"/><path d="M10 2v3h3"/>
          </svg>
          <span className="path"><b>{yamlName}</b></span>
        </div>
        <button className="btn ghost sm" onClick={() => {
          const el = document.createElement("input");
          el.type = "file"; el.accept = ".yaml,.yml";
          el.onchange = () => el.files[0] && handleLoadFile(el.files[0]);
          el.click();
        }}>Open…</button>
        <button className="btn ghost sm" onClick={() => setEditorOpen(true)}
                disabled={!yamlText}>Edit YAML…</button>
        <button className="btn ghost sm" onClick={handleRun}
                disabled={!yamlText || loading}>Reload</button>
        <div className="topbar-spacer"/>
        <div style={{display:"flex",alignItems:"center",gap:6, color:"var(--ink-3)", fontSize:11, fontFamily:"var(--ff-mono)"}}>
          <span className="dot" style={{
            width:6, height:6, borderRadius:"50%",
            background: loading ? "var(--warn)" : (error ? "var(--fail)" : "var(--ok)"),
            display:"inline-block",
          }}/>
          {loading ? "solving…" : (error ? "error" : "idle")}
        </div>
        <div className="topbar-sep"/>
        <button className="btn ghost icon" title="Theme"
                onClick={()=>setTweak("dark", !tweaks.dark)}>
          {tweaks.dark
            ? <svg className="ic" viewBox="0 0 16 16" fill="currentColor"><path d="M6 2a6 6 0 1 0 8 8A7 7 0 0 1 6 2z"/></svg>
            : <svg className="ic" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.4"><circle cx="8" cy="8" r="3"/><path d="M8 1v2M8 13v2M1 8h2M13 8h2M3 3l1.5 1.5M11.5 11.5L13 13M3 13l1.5-1.5M11.5 4.5L13 3"/></svg>}
        </button>
        <button className="btn primary" onClick={handleRun}
                disabled={!yamlText || loading}>
          {loading ? "Solving…" : "▸ Run"}
        </button>
      </div>

      {/* MAIN */}
      <div className="main">
        <LeftRail selected={selected} onSelect={setSelected} />

        <section className="center">
          <div className="viewport">
            <div className="vp-tabs">
              {tabs.map(t => (
                <button key={t.id} className={"vp-tab" + (tab===t.id ? " active":"")}
                        onClick={()=>setTab(t.id)} disabled={isEmpty}>{t.label}</button>
              ))}
              <div className="vp-spacer"/>
              <div className="vp-tools">
                <button className="btn ghost sm" disabled>PNG</button>
                <button className="btn ghost sm" disabled>CSV</button>
                <button className="btn ghost sm" disabled>JSON</button>
              </div>
            </div>

            <div className="vp-stage" ref={stageRef}>
              {isEmpty ? (
                <EmptyState
                  onLoadFile={handleLoadFile}
                  onLoadExample={handleLoadExample}
                  loading={loading}
                  error={error}
                />
              ) : (
                <>
                  <div className="plot-caption">{tabs.find(t=>t.id===tab).cap}</div>
                  <div className="plot-subcaption">
                    {tab === "mxmy"    && `N = ${nSlice} kN   ·   section ${D.SECTION.B}×${D.SECTION.H}`}
                    {tab === "nm"      && `Uniaxial about X   ·   section ${D.SECTION.B}×${D.SECTION.H}`}
                    {tab === "surface" && `convex hull   ·   ${D.MATERIALS.map(m=>m.id).join(" / ")}`}
                    {tab === "mchi"    && `per N-level   ·   ${D.MATERIALS.map(m=>m.id).join(" / ")}`}
                    {tab === "polar"   && "ductility ratio μ_χ"}
                    {tab === "section" && `${D.SECTION.n_fibers_x}×${D.SECTION.n_fibers_y} fibers   ·   ${D.SECTION.rebars.length} rebars`}
                  </div>

                  {tab === "surface" && <Surface3D width={stageSize.w} height={stageSize.h} />}
                  {tab === "mxmy" && (
                    <>
                      <MxMyPlot width={stageSize.w} height={stageSize.h - 56}
                                N_kN={nSlice} demands={demandPoints} activeName={activeName}/>
                      <div className="n-slider">
                        <span className="label">N slice</span>
                        <input type="range" min={-4500} max={1200} step={50}
                               value={nSlice}
                               onChange={(e)=>setNSlice(parseInt(e.target.value))}/>
                        <span className="value">{nSlice} kN</span>
                      </div>
                    </>
                  )}
                  {tab === "nm"      && <NMPlot width={stageSize.w} height={stageSize.h}
                                                demands={demandPoints} activeName={activeName}/>}
                  {tab === "mchi"    && <MchiPlot width={stageSize.w} height={stageSize.h} activeN={nSlice}/>}
                  {tab === "polar"   && <PolarDuctility width={stageSize.w} height={stageSize.h}/>}
                  {tab === "section" && <SectionPreview width={stageSize.w} height={stageSize.h}/>}

                  {(tab === "mxmy" || tab === "nm" || tab === "surface") && (
                    <div className="legend">
                      <div className="row">
                        <span className="sw" style={{background:"var(--accent)"}}/> Resistance domain
                      </div>
                      <div className="row">
                        <span className="sw" style={{background:"var(--ink)", borderRadius:"50%"}}/> Demand
                      </div>
                      <div className="row">
                        <span className="sw" style={{background:"var(--fail)", borderRadius:"50%"}}/> Selected · η-ray
                      </div>
                    </div>
                  )}
                </>
              )}
            </div>
          </div>

          <div className="results">
            <div className="results-tabs">
              {resTabs.map(t => (
                <button key={t.id}
                        className={"results-tab" + (resTab===t.id ? " active":"")}
                        onClick={()=>setResTab(t.id)}>
                  {t.label}<span className="count">{t.count}</span>
                </button>
              ))}
              <div className="vp-spacer"/>
              <div style={{display:"flex", gap:6, padding: "0 6px"}}>
                <button className="btn ghost sm" disabled>Filter</button>
                <button className="btn ghost sm" disabled>Export summary</button>
              </div>
            </div>
            <div className="results-body">
              {isEmpty ? (
                <div style={{
                  padding: 32, textAlign: "center", color: "var(--ink-3)",
                  fontSize: 12, fontFamily: "var(--ff-mono)",
                }}>
                  Load a YAML to see results.
                </div>
              ) : resTab === "materials" ? (
                <table className="dt">
                  <thead><tr>
                    <th style={{width:24}}></th><th>ID</th><th>Type</th><th>Class</th>
                    <th>Design strength</th><th>Modulus</th><th>Ultimate</th>
                  </tr></thead>
                  <tbody>
                    {D.MATERIALS.map(m => (
                      <tr key={m.id}>
                        <td><span className="dot" style={{background:"var(--series-2)"}}/></td>
                        <td className="mono" style={{fontWeight:500}}>{m.id}</td>
                        <td style={{color:"var(--ink-3)"}}>{m.kind}</td>
                        <td className="mono">{m.cls}</td>
                        <td className="mono">{m.fcd || m.fyd || "—"}</td>
                        <td className="mono">{m.Ecm || m.Es || "—"}</td>
                        <td className="mono">{m.eps_su || "—"}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              ) : (
                <ResultsTable rows={resRows} selected={selected}
                              onSelect={setSelected} flags={flags}/>
              )}
            </div>
          </div>
        </section>

        <RightDrawer flags={flags} setFlag={setFlag}
                     nSlice={nSlice} setNSlice={setNSlice}
                     onRun={handleRun} />
      </div>

      {/* STATUSBAR */}
      <div className="statusbar">
        <span>
          <span className="dot" style={{background: loading ? "var(--warn)" : "var(--ok)"}}/>
          gensec {loading ? "solving" : "ready"}
        </span>
        <span>{yamlName}</span>
        <span>selected · <b style={{color:"var(--ink-2)"}}>{selected || "—"}</b></span>
        <div className="spacer"/>
        <span>{D.MATERIALS.length} materials · {D.DEMANDS.length} demands · {D.COMBINATIONS.length} combinations · {D.ENVELOPES.length} envelopes</span>
        {D._meta && D._meta.elapsed_ms !== undefined &&
          <span>solved in {D._meta.elapsed_ms} ms{D._meta.cached ? " (cached)" : ""}</span>}
      </div>

      <TweaksPanel open={tweaksOpen} onClose={()=>setTweaksOpen(false)}
                   tweaks={tweaks} setTweak={setTweak} />

      <YamlEditor
        open={editorOpen}
        initialText={yamlText}
        initialName={yamlName}
        onClose={() => setEditorOpen(false)}
        onSubmit={(text, name) => {
          setEditorOpen(false);
          runAnalyze(text, name);
        }}
      />

      <ErrorBanner message={!isEmpty && error} onDismiss={() => setError(null)} />
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App/>);
