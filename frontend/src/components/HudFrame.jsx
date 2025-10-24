const HudFrame = ({ id, title, subtitle, actions, children, className = '' }) => (
  <div
    id={id}
    className={`flex flex-col rounded-xl border border-orange-200 bg-white shadow-sm ${className}`}
  >
    {(title || subtitle || actions) && (
      <div className="flex items-start justify-between border-b border-orange-200 px-4 py-3">
        <div>
          {title && <h3 className="text-sm font-semibold text-orange-700">{title}</h3>}
          {subtitle && <p className="text-xs text-orange-600/70">{subtitle}</p>}
        </div>
        {actions && <div className="flex items-center gap-2 text-sm text-orange-600/80">{actions}</div>}
      </div>
    )}
    <div className="flex-1 overflow-hidden px-4 py-3">
      {children}
    </div>
  </div>
);

export default HudFrame;
