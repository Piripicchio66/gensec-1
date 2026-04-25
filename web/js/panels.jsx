/* ==========================================================
 * GenSec GUI — side panels: left rail tree, right drawer, results table
 * ========================================================== */

function LeftRail({ selected, onSelect }) {
  const D = window.GS_DATA;
  const Section = ({ title, extra, children }) => (
    <div className="tree-section">
      <div className="overline">
        <span>{title}</span>
        <span>{extra}</span>
      </div>
      {children}
    </div>
  );
  const Row = ({ kind, name, aux, id }) => (
    <div className={"tree-row " + kind + (selected === id ? " active" : "")}
         onClick={() => onSelect(id)}>
      <span className="bullet" />
      <span className="name">{name}</span>
      {aux && <span className="aux">{aux}</span>}
    </div>
  );
  return (
    <aside className="rail">
      <div className="panel-head">
        <h2>Input · YAML</h2>
        <span className="meta">example_v2_1</span>
      </div>
      <div className="panel-body">
        <Section title="Materials" extra={D.MATERIALS.length}>
          {D.MATERIALS.map(m => (
            <Row key={m.id} kind="material" id={"mat:"+m.id}
                 name={m.id} aux={m.cls} />
          ))}
        </Section>
        <Section title="Section" extra={`${D.SECTION.B}×${D.SECTION.H}`}>
          <Row kind="material" id="section" name="geometry"
               aux={`${D.SECTION.n_fibers_x}×${D.SECTION.n_fibers_y}`} />
          <Row kind="material" id="rebars"  name="rebars"
               aux={D.SECTION.rebars.length + " bars"} />
        </Section>
        <Section title="Demands" extra={D.DEMANDS.length}>
          {D.DEMANDS.map(d => (
            <Row key={d.name} kind="demand" id={"dem:"+d.name}
                 name={d.name} aux={`N ${d.N}`} />
          ))}
        </Section>
        <Section title="Combinations" extra={D.COMBINATIONS.length}>
          {D.COMBINATIONS.map(c => (
            <Row key={c.name} kind="combo" id={"cmb:"+c.name}
                 name={c.name} aux={c.staged ? "staged" : "simple"} />
          ))}
        </Section>
        <Section title="Envelopes" extra={D.ENVELOPES.length}>
          {D.ENVELOPES.map(e => (
            <Row key={e.name} kind="envelope" id={"env:"+e.name}
                 name={e.name} aux={e.members.length + " mbrs"} />
          ))}
        </Section>
      </div>
    </aside>
  );
}

function Toggle({ on, onChange }) {
  return <div className={"sw-toggle" + (on ? " on" : "")} onClick={() => onChange(!on)} />;
}

function RightDrawer({ flags, setFlag, nSlice, setNSlice, onRun }) {
  return (
    <aside className="drawer">
      <div className="panel-head">
        <h2>Output flags</h2>
        <span className="meta">v2.1</span>
      </div>
      <div className="panel-body">
        <div className="drawer-section">
          <h3>Utilization ratios</h3>
          {[
            ["eta_3D", "3D hull ray from origin"],
            ["eta_2D", "2D ray on Mx–My contour"],
            ["eta_path", "3D staged ray"],
            ["eta_path_2D", "2D staged ray"],
          ].map(([k, desc]) => (
            <div key={k} className="flag-row">
              <span className="lbl"><code>{k}</code></span>
              <Toggle on={flags[k]} onChange={(v) => setFlag(k, v)} />
            </div>
          ))}
          <div className="field">
            <label>delta_N_tol</label>
            <input type="number" value={flags.delta_N_tol} step="0.005"
                   onChange={(e)=>setFlag("delta_N_tol", parseFloat(e.target.value))}/>
          </div>
        </div>

        <div className="drawer-section">
          <h3>Domain generation</h3>
          {[
            ["generate_mx_my", "Mx–My contours"],
            ["generate_3d_surface", "3D resistance surface"],
            ["generate_moment_curvature", "M–χ diagrams"],
            ["generate_polar_ductility", "Polar ductility"],
            ["generate_3d_moment_curvature", "3D M–χ–N surface"],
          ].map(([k]) => (
            <div key={k} className="flag-row">
              <span className="lbl"><code>{k}</code></span>
              <Toggle on={flags[k]} onChange={(v) => setFlag(k, v)} />
            </div>
          ))}
          <div className="field">
            <label>n_angles_mx_my</label>
            <input type="number" value={flags.n_angles_mx_my} step="12"
                   onChange={(e)=>setFlag("n_angles_mx_my", parseInt(e.target.value))}/>
          </div>
        </div>

        <div className="drawer-section">
          <h3>Viewport · N slice</h3>
          <div className="field">
            <label>N  [kN]</label>
            <span className="mono" style={{minWidth:70, textAlign:"right"}}>{nSlice}</span>
          </div>
          <input type="range" min={-4500} max={1200} step={50}
                 value={nSlice} onChange={(e)=>setNSlice(parseInt(e.target.value))}
                 style={{width: "100%", accentColor: "var(--accent)"}} />
        </div>

        <div className="drawer-section">
          <h3>Solver</h3>
          <div className="field"><label>mesh_method</label>
            <select defaultValue="grid">
              <option>grid</option><option>triangle</option>
            </select></div>
          <div className="field"><label>n_fibers_y</label>
            <input type="number" defaultValue={40}/></div>
          <div className="field"><label>n_fibers_x</label>
            <input type="number" defaultValue={20}/></div>
          <button className="btn primary" style={{width:"100%", marginTop:8, justifyContent:"center"}}
                  onClick={onRun}>
            ▸ Run gensec
          </button>
        </div>
      </div>
    </aside>
  );
}

function EtaPill({ v }) {
  if (v === null || v === undefined) return <span style={{color:"var(--ink-3)"}}>—</span>;
  const cls = v >= 1 ? "fail" : v >= 0.85 ? "warn" : "ok";
  return <span className={"pill " + cls}>η = {v.toFixed(2)}</span>;
}

function ResultsTable({ rows, selected, onSelect, flags }) {
  return (
    <table className="dt">
      <thead>
        <tr>
          <th style={{width: 24}}></th>
          <th style={{minWidth: 180}}>Name</th>
          <th>Type</th>
          <th className="num">N [kN]</th>
          <th className="num">Mx [kN·m]</th>
          <th className="num">My [kN·m]</th>
          {flags.eta_3D     && <th className="num">η₃D</th>}
          {flags.eta_2D     && <th className="num">η₂D</th>}
          {flags.eta_path   && <th className="num">η<sub>path</sub></th>}
          {flags.eta_path_2D&& <th className="num">η<sub>path,2D</sub></th>}
          <th>Status</th>
        </tr>
      </thead>
      <tbody>
        {rows.map((r) => {
          const id = (r.kind === "demand" ? "dem:" : r.kind === "combo" ? "cmb:" : "env:") + r.name;
          const active = selected === id;
          const dot = r.kind === "demand" ? "var(--series-1)"
                    : r.kind === "combo"  ? "var(--series-5)"
                                          : "var(--series-3)";
          return (
            <tr key={id} className={active ? "active" : ""} onClick={()=>onSelect(id)}>
              <td><span className="dot" style={{background: dot}}/></td>
              <td className="mono" style={{fontWeight: 500}}>{r.name}</td>
              <td style={{color:"var(--ink-3)", fontSize: 11}}>{r.kind}</td>
              <td className="num">{typeof r.N === "number" ? r.N : r.N}</td>
              <td className="num">{typeof r.Mx === "number" ? r.Mx : r.Mx}</td>
              <td className="num">{typeof r.My === "number" ? r.My : r.My}</td>
              {flags.eta_3D     && <td className="num"><EtaPill v={r.eta3D}/></td>}
              {flags.eta_2D     && <td className="num"><EtaPill v={r.eta2D}/></td>}
              {flags.eta_path   && <td className="num"><EtaPill v={r.etaPath}/></td>}
              {flags.eta_path_2D&& <td className="num"><EtaPill v={r.etaPath2D}/></td>}
              <td>
                {r.status === "ok"   && <span className="pill ok">✓ verified</span>}
                {r.status === "warn" && <span className="pill warn">⚠ close</span>}
                {r.status === "fail" && <span className="pill fail">✗ exceeds</span>}
              </td>
            </tr>
          );
        })}
      </tbody>
    </table>
  );
}

Object.assign(window, { LeftRail, RightDrawer, ResultsTable, Toggle });
