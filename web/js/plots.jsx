/* ==========================================================
 * GenSec GUI — SVG plots
 * All plots are static-friendly but accept {active} highlighting.
 * ========================================================== */

// ----- utilities -----
function axisTicks(min, max, count = 6) {
  const step = (max - min) / (count - 1);
  const ticks = [];
  for (let i = 0; i < count; i++) ticks.push(min + i * step);
  return ticks;
}
function fmt(n, d = 0) {
  if (n === null || n === undefined || isNaN(n)) return "—";
  if (Math.abs(n) >= 1000) return n.toFixed(d).replace(/\B(?=(\d{3})+(?!\d))/g, " ");
  return n.toFixed(d);
}

function PlotFrame({ width, height, pad, xDom, yDom, xLabel, yLabel, children, grid=true }) {
  const [x0, x1] = xDom, [y0, y1] = yDom;
  const sx = (v) => pad.l + ((v - x0) / (x1 - x0)) * (width - pad.l - pad.r);
  const sy = (v) => (height - pad.b) - ((v - y0) / (y1 - y0)) * (height - pad.t - pad.b);
  const xTicks = axisTicks(x0, x1, 7);
  const yTicks = axisTicks(y0, y1, 6);
  return (
    <svg width={width} height={height} style={{ display: "block" }}>
      {/* grid */}
      {grid && xTicks.map((t, i) => (
        <line key={"gx"+i} x1={sx(t)} x2={sx(t)} y1={pad.t} y2={height - pad.b} stroke="var(--rule)" strokeDasharray="1 3" />
      ))}
      {grid && yTicks.map((t, i) => (
        <line key={"gy"+i} x1={pad.l} x2={width - pad.r} y1={sy(t)} y2={sy(t)} stroke="var(--rule)" strokeDasharray="1 3" />
      ))}
      {/* axes */}
      <line x1={pad.l} x2={width - pad.r} y1={height - pad.b} y2={height - pad.b} stroke="var(--rule-2)" />
      <line x1={pad.l} x2={pad.l} y1={pad.t} y2={height - pad.b} stroke="var(--rule-2)" />
      {/* origin */}
      {x0 < 0 && x1 > 0 && (
        <line x1={sx(0)} x2={sx(0)} y1={pad.t} y2={height - pad.b} stroke="var(--ink-3)" strokeWidth="0.6" />
      )}
      {y0 < 0 && y1 > 0 && (
        <line x1={pad.l} x2={width - pad.r} y1={sy(0)} y2={sy(0)} stroke="var(--ink-3)" strokeWidth="0.6" />
      )}
      {/* tick labels */}
      {xTicks.map((t, i) => (
        <text key={"tx"+i} x={sx(t)} y={height - pad.b + 14} textAnchor="middle"
              fontFamily="var(--ff-mono)" fontSize="10" fill="var(--ink-3)">{fmt(t)}</text>
      ))}
      {yTicks.map((t, i) => (
        <text key={"ty"+i} x={pad.l - 8} y={sy(t) + 3} textAnchor="end"
              fontFamily="var(--ff-mono)" fontSize="10" fill="var(--ink-3)">{fmt(t)}</text>
      ))}
      {/* axis labels */}
      <text x={width - pad.r} y={height - 6} textAnchor="end"
            fontFamily="var(--ff-serif)" fontStyle="italic" fontSize="11" fill="var(--ink-2)">{xLabel}</text>
      <text x={pad.l - 4} y={pad.t - 6} textAnchor="start"
            fontFamily="var(--ff-serif)" fontStyle="italic" fontSize="11" fill="var(--ink-2)">{yLabel}</text>
      {children(sx, sy)}
    </svg>
  );
}

/* ========== Mx-My contour plot at fixed N ========== */
function MxMyPlot({ width, height, N_kN, demands, activeName }) {
  const contour = window.GS_DATA.makeMxMyContour(N_kN, 144);
  const pad = { l: 52, r: 20, t: 38, b: 32 };

  // auto-extents with padding
  const xs = contour.map(p => p[0]); const ys = contour.map(p => p[1]);
  const xmax = Math.max(420, Math.max(...xs.map(Math.abs)) * 1.15);
  const ymax = Math.max(280, Math.max(...ys.map(Math.abs)) * 1.15);

  return (
    <PlotFrame width={width} height={height} pad={pad}
               xDom={[-xmax, xmax]} yDom={[-ymax, ymax]}
               xLabel="Mx  [kN·m]" yLabel="My  [kN·m]">
      {(sx, sy) => (
        <>
          {/* filled contour */}
          <path d={"M " + contour.map(p => sx(p[0]) + " " + sy(p[1])).join(" L ") + " Z"}
                fill="var(--accent-soft)" fillOpacity="0.35"
                stroke="var(--accent)" strokeWidth="1.6" />
          {/* demand points */}
          {demands.map((d, i) => {
            const active = d.name === activeName;
            const cx = sx(d.Mx), cy = sy(d.My);
            // ray from origin to contour through demand direction, drawn if active
            return (
              <g key={i}>
                {active && (
                  <line x1={sx(0)} y1={sy(0)} x2={cx} y2={cy}
                        stroke="var(--fail)" strokeWidth="1.4" strokeDasharray="3 2" />
                )}
                <circle cx={cx} cy={cy} r={active ? 5 : 3.2}
                        fill={active ? "var(--fail)" : "var(--ink)"}
                        stroke="var(--paper)" strokeWidth="1.5" />
                <text x={cx + 7} y={cy - 6} fontFamily="var(--ff-mono)" fontSize="10"
                      fill={active ? "var(--fail)" : "var(--ink-2)"} fontWeight={active ? 600 : 400}>
                  {d.name}
                </text>
              </g>
            );
          })}
        </>
      )}
    </PlotFrame>
  );
}

/* ========== N-M (uniaxial) interaction ========== */
function NMPlot({ width, height, demands, activeName }) {
  const pts = window.GS_DATA.makeNM();
  const pad = { l: 60, r: 20, t: 38, b: 32 };
  return (
    <PlotFrame width={width} height={height} pad={pad}
               xDom={[-5000, 1500]} yDom={[-460, 460]}
               xLabel="N  [kN]" yLabel="Mx  [kN·m]">
      {(sx, sy) => (
        <>
          <path d={"M " + pts.map(p => sx(p[0]) + " " + sy(p[1])).join(" L ") + " Z"}
                fill="var(--accent-soft)" fillOpacity="0.35"
                stroke="var(--accent)" strokeWidth="1.6" />
          {demands.map((d, i) => {
            const active = d.name === activeName;
            const cx = sx(d.N), cy = sy(d.Mx);
            return (
              <g key={i}>
                <circle cx={cx} cy={cy} r={active ? 5 : 3.2}
                        fill={active ? "var(--fail)" : "var(--ink)"}
                        stroke="var(--paper)" strokeWidth="1.5" />
                <text x={cx + 7} y={cy - 6} fontFamily="var(--ff-mono)" fontSize="10"
                      fill={active ? "var(--fail)" : "var(--ink-2)"} fontWeight={active ? 600 : 400}>
                  {d.name}
                </text>
              </g>
            );
          })}
        </>
      )}
    </PlotFrame>
  );
}

/* ========== 3D Resistance Surface (isometric stack of slices) ========== */
function Surface3D({ width, height }) {
  const slices = window.GS_DATA.makeSurfaceSlices(11);
  // map (Mx, My, N) → 2D via axonometric projection
  const ang = 28 * Math.PI / 180;
  const cx = width / 2;
  const cy = height / 2 + 20;
  const sMx = 0.35;        // Mx scale
  const sMy = 0.35;        // My scale
  const sN  = 0.055;       // N scale (vertical)
  const proj = (Mx, My, N) => {
    // classic isometric-ish: Mx along cos(30), -sin(30); My along cos(-30), -sin(-30)
    const x = cx + Mx * sMx * Math.cos(ang) + My * sMy * Math.cos(Math.PI - ang);
    const y = cy + Mx * sMx * Math.sin(ang) + My * sMy * Math.sin(Math.PI - ang) - N * sN;
    return [x, y];
  };

  // Build stacked contours as polygons, colored by N level
  return (
    <svg width={width} height={height} style={{ display: "block" }}>
      {/* vertical N axis guide */}
      <line x1={cx} y1={cy - (-4500) * sN * -1} x2={cx} y2={cy - 1200 * sN * -1}
            stroke="var(--rule-2)" strokeDasharray="2 3" />
      {/* floor ellipses at a few N-levels as faint shadow */}
      {slices.map((s, i) => {
        const pts = s.contour.map(p => proj(p[0], p[1], s.N));
        const d = "M " + pts.map(p => p[0].toFixed(1) + " " + p[1].toFixed(1)).join(" L ") + " Z";
        const t = i / (slices.length - 1);
        const hueL = 30 + t * 50;
        return (
          <path key={i} d={d}
                fill={`oklch(${hueL}% 0.12 225 / 0.18)`}
                stroke="var(--accent)" strokeWidth="1.1" strokeOpacity={0.55 + t * 0.3} />
        );
      })}
      {/* axes labels */}
      <g fontFamily="var(--ff-serif)" fontStyle="italic" fontSize="11" fill="var(--ink-2)">
        {(() => {
          const [ax, ay] = proj(420, 0, 0);
          const [bx, by] = proj(0, 280, 0);
          const [nx, ny] = proj(0, 0, 1500);
          return (
            <>
              <line x1={proj(0,0,0)[0]} y1={proj(0,0,0)[1]} x2={ax} y2={ay} stroke="var(--ink-3)" markerEnd=""/>
              <text x={ax + 6} y={ay + 4}>Mx</text>
              <line x1={proj(0,0,0)[0]} y1={proj(0,0,0)[1]} x2={bx} y2={by} stroke="var(--ink-3)"/>
              <text x={bx - 22} y={by + 4}>My</text>
              <line x1={cx} y1={cy - (-4500) * -sN} x2={cx} y2={cy - 1500 * -sN} stroke="var(--ink-3)"/>
              <text x={cx + 6} y={cy - 1500 * -sN - 6}>N</text>
            </>
          );
        })()}
      </g>
      {/* N tick labels on central axis */}
      {[-4000, -3000, -2000, -1000, 0, 1000].map((N, i) => {
        const [x, y] = proj(0, 0, N);
        return (
          <g key={i}>
            <line x1={x - 3} y1={y} x2={x + 3} y2={y} stroke="var(--ink-3)"/>
            <text x={x - 6} y={y + 3} textAnchor="end"
                  fontFamily="var(--ff-mono)" fontSize="9.5" fill="var(--ink-3)">{N}</text>
          </g>
        );
      })}
    </svg>
  );
}

/* ========== M-χ curves ========== */
function MchiPlot({ width, height, activeN }) {
  const series = window.GS_DATA.makeMchi();
  const pad = { l: 56, r: 20, t: 38, b: 32 };
  // scale curvature to 1e-6 /mm for readability
  return (
    <PlotFrame width={width} height={height} pad={pad}
               xDom={[0, 145]} yDom={[0, 440]}
               xLabel="χ  ×10⁻⁶  [1/mm]" yLabel="M  [kN·m]">
      {(sx, sy) => (
        <>
          {series.map((s, i) => {
            const active = activeN === s.N;
            const color = ["var(--series-1)","var(--series-2)","var(--series-3)","var(--series-4)","var(--series-5)"][i % 5];
            const d = s.pts.map((p, j) => (j === 0 ? "M " : "L ") + sx(p[0]*1e6) + " " + sy(p[1])).join(" ");
            return <path key={i} d={d} fill="none" stroke={color}
                         strokeWidth={active ? 2.2 : 1.3}
                         strokeOpacity={active ? 1 : 0.7} />;
          })}
        </>
      )}
    </PlotFrame>
  );
}

/* ========== Section preview with fiber mesh + rebars ========== */
function SectionPreview({ width, height }) {
  const { B, H, rebars, n_fibers_x, n_fibers_y } = window.GS_DATA.SECTION;
  const pad = 32;
  const sc = Math.min((width - pad*2) / B, (height - pad*2) / H);
  const w = B * sc, h = H * sc;
  const x0 = (width - w) / 2, y0 = (height - h) / 2;

  const dx = w / n_fibers_x, dy = h / n_fibers_y;
  const lines = [];
  for (let i = 1; i < n_fibers_x; i++) lines.push(<line key={"vx"+i} x1={x0+i*dx} x2={x0+i*dx} y1={y0} y2={y0+h} stroke="var(--rule)" strokeWidth="0.5"/>);
  for (let j = 1; j < n_fibers_y; j++) lines.push(<line key={"vy"+j} x1={x0} x2={x0+w} y1={y0+j*dy} y2={y0+j*dy} stroke="var(--rule)" strokeWidth="0.5"/>);

  return (
    <svg width={width} height={height} style={{ display: "block" }}>
      {/* outline */}
      <rect x={x0} y={y0} width={w} height={h}
            fill="oklch(88% 0.015 85 / 0.5)" stroke="var(--ink-2)" strokeWidth="1.2" />
      {lines}
      {/* rebars (y measured from bottom — invert to SVG) */}
      {rebars.map((r, i) => {
        const cx = x0 + (r.x / B) * w;
        const cy = y0 + h - (r.y / H) * h;
        const rad = (r.diameter / 2) * sc * 1.6;
        return (
          <g key={i}>
            <circle cx={cx} cy={cy} r={rad} fill="var(--ink)" />
            <circle cx={cx} cy={cy} r={rad*0.5} fill="var(--paper)" opacity="0.15" />
          </g>
        );
      })}
      {/* dims */}
      <g fontFamily="var(--ff-mono)" fontSize="10" fill="var(--ink-3)">
        <text x={x0 + w/2} y={y0 + h + 16} textAnchor="middle">B = {B} mm</text>
        <text x={x0 - 10} y={y0 + h/2} textAnchor="end" transform={`rotate(-90 ${x0-10} ${y0+h/2})`}>H = {H} mm</text>
      </g>
      <g fontFamily="var(--ff-serif)" fontStyle="italic" fontSize="11" fill="var(--ink-2)">
        <text x={width - 14} y={20} textAnchor="end">fiber mesh {n_fibers_x} × {n_fibers_y}</text>
      </g>
    </svg>
  );
}

/* ========== Polar ductility (simple rose) ========== */
function PolarDuctility({ width, height }) {
  const cx = width/2, cy = height/2 + 6;
  const R = Math.min(width, height) / 2 - 40;
  const n = 36;
  const pts = [];
  for (let i = 0; i <= n; i++) {
    const t = (i/n) * Math.PI * 2;
    const r = R * (0.55 + 0.32 * Math.sin(2*t + 0.3) + 0.08 * Math.cos(4*t));
    pts.push([cx + Math.cos(t)*r, cy + Math.sin(t)*r]);
  }
  return (
    <svg width={width} height={height} style={{display:"block"}}>
      {[0.25, 0.5, 0.75, 1].map((f, i) => (
        <circle key={i} cx={cx} cy={cy} r={R*f} fill="none" stroke="var(--rule)" strokeDasharray="1 3" />
      ))}
      {[0, 45, 90, 135, 180, 225, 270, 315].map((a, i) => {
        const rad = a * Math.PI/180;
        return <line key={i} x1={cx} y1={cy} x2={cx+Math.cos(rad)*R} y2={cy+Math.sin(rad)*R}
                     stroke="var(--rule)" strokeDasharray="1 3" />;
      })}
      <path d={"M " + pts.map(p => p[0].toFixed(1)+" "+p[1].toFixed(1)).join(" L ") + " Z"}
            fill="var(--accent-soft)" fillOpacity="0.5"
            stroke="var(--accent)" strokeWidth="1.6"/>
      <text x={cx + R + 6} y={cy + 4} fontFamily="var(--ff-mono)" fontSize="10" fill="var(--ink-3)">Mx</text>
      <text x={cx - 4} y={cy - R - 6} fontFamily="var(--ff-mono)" fontSize="10" fill="var(--ink-3)">My</text>
    </svg>
  );
}

Object.assign(window, { MxMyPlot, NMPlot, Surface3D, MchiPlot, SectionPreview, PolarDuctility });
