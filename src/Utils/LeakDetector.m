function report = LeakDetector(details, q_in, varargin)
% LeakDetector  Simple per-run leak detection heuristic
% Usage:
%   report = LeakDetector(details, q_in)
%   report = LeakDetector(details, q_in, 'tol_dp',1000, 'imbalance_frac',0.01)

% parse options
p = inputParser;
addParameter(p,'tol_dp',1e3,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'imbalance_frac',0.01,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'min_severity_threshold',0.05,@(x)isnumeric(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;

% If solver residual is huge, suppress leak flag (avoid false positives)
if isfield(details,'res_norm') && ~isempty(details.res_norm) && isfinite(details.res_norm) && details.res_norm > 1e3
  report = struct('leak_flag',false,'severity',0,'suspected_segment',NaN,'dp_seg',getfield_or_empty(details,'dp_seg'),'p_seg',getfield_or_empty(details,'p_seg'),'imbalance_frac',NaN,'timestamp',datetime('now'));
  return;
end

% ensure dp_seg and p_seg exist
if ~isfield(details,'dp_seg') || ~isfield(details,'p_seg')
  if isfield(details,'dp_total')
    Nseg = 20;
    dp_seg = repmat(details.dp_total / Nseg, 1, Nseg);
    if isfield(details,'BHP') && isfield(details,'p_surface_pred')
      p_seg = linspace(details.BHP, details.p_surface_pred, Nseg+1);
    else
      p_seg = nan(1,Nseg+1);
    end
  else
    report = struct('leak_flag',false,'severity',0,'suspected_segment',NaN,'dp_seg',[],'p_seg',[],'imbalance_frac',NaN,'timestamp',datetime('now'));
    return;
  end
else
  dp_seg = details.dp_seg;
  p_seg = details.p_seg;
end

% basic anomaly detection: per-segment dp compared to tol_dp
dp_abs = abs(dp_seg);
if isempty(dp_abs)
  report = struct('leak_flag',false,'severity',0,'suspected_segment',NaN,'dp_seg',dp_seg,'p_seg',p_seg,'imbalance_frac',NaN,'timestamp',datetime('now'));
  return;
end

% normalized anomaly magnitude
norm_anom = dp_abs / max(opts.tol_dp,1);
norm_anom(~isfinite(norm_anom)) = 0;

% contrast rule: require anomaly to be significantly larger than local/median dp
median_dp = median(dp_abs + eps);
contrast_mask = dp_abs > 5 * median_dp; % 5x median by default

% severity raw: average normalized anomaly weighted by contrast mask
if any(contrast_mask)
  severity_raw = sum(norm_anom .* double(contrast_mask)) / sum(double(contrast_mask));
else
  severity_raw = 0;
end

% mass imbalance proxy: placeholder (kept for future extension)
imbalance_frac = opts.imbalance_frac;
if isfield(details,'imbalance_frac') && ~isempty(details.imbalance_frac) && ~isnan(details.imbalance_frac)
  imbalance_frac = details.imbalance_frac;
end
if ~isfinite(imbalance_frac), imbalance_frac = 0; end
imbalance_frac = max(0, min(1, imbalance_frac));

% combine metrics into severity 0..1 (80% anomaly, 20% imbalance)
severity = min(1, 0.8 * severity_raw + 0.2 * imbalance_frac);

% apply minimum threshold
leak_flag = severity >= opts.min_severity_threshold;

% localize: pick segment with largest normalized anomaly
[~, suspected_idx] = max(norm_anom);

report.leak_flag = leak_flag;
report.severity = severity;
report.suspected_segment = suspected_idx;
report.dp_seg = dp_seg;
report.p_seg = p_seg;
report.imbalance_frac = imbalance_frac;
report.timestamp = datetime('now');
end

% small helper to safely get a field or return empty
function v = getfield_or_empty(s, name)
if isstruct(s) && isfield(s,name)
  v = s.(name);
else
  v = [];
end
end
