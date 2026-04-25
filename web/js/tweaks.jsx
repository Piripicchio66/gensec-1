/* Tweaks floating panel */
function TweaksPanel({ open, onClose, tweaks, setTweak }) {
  if (!open) return null;
  const huePct = ((tweaks.accent_h % 360) / 360) * 100;
  return (
    <div className="tweaks">
      <div className="tweaks-head">
        <h3>Tweaks</h3>
        <button className="btn ghost sm" onClick={onClose}>×</button>
      </div>
      <div className="tweaks-body">
        <div className="overline" style={{marginBottom:4}}>Accent hue</div>
        <div className="hue-strip" onClick={(e)=>{
          const r = e.currentTarget.getBoundingClientRect();
          const p = (e.clientX - r.left) / r.width;
          setTweak("accent_h", Math.round(p * 360));
        }}>
          <div className="knob" style={{left: `${huePct}%`}} />
        </div>

        <div className="flag-row">
          <span className="lbl">Dark theme</span>
          <Toggle on={tweaks.dark} onChange={(v)=>setTweak("dark", v)} />
        </div>
        <div className="flag-row">
          <span className="lbl">Serif headings</span>
          <Toggle on={tweaks.serif} onChange={(v)=>setTweak("serif", v)} />
        </div>
        <div className="flag-row">
          <span className="lbl">Compact density</span>
          <Toggle on={tweaks.compact} onChange={(v)=>setTweak("compact", v)} />
        </div>
        <div className="flag-row">
          <span className="lbl">Show grid on plots</span>
          <Toggle on={tweaks.grid} onChange={(v)=>setTweak("grid", v)} />
        </div>

        <div className="field" style={{marginTop:8}}>
          <label>Left rail width</label>
          <input type="number" value={tweaks.rail_w} min={200} max={360} step={10}
                 onChange={(e)=>setTweak("rail_w", parseInt(e.target.value))}/>
        </div>
      </div>
    </div>
  );
}
window.TweaksPanel = TweaksPanel;
