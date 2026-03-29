function validate_leak_detector(outFile)
% validate_leak_detector  Synthetic injection tests across depths and magnitudes
% Usage: validate_leak_detector() or validate_leak_detector('results/validation_results.csv')

if nargin < 1 || isempty(outFile)
  outFile = fullfile(pwd,'results','validation_results.csv');
end

% determine project root robustly
thisFile = mfilename('fullpath');
if isempty(thisFile)
  proj = pwd;
else
  proj = fileparts(fileparts(thisFile));
end
addpath(genpath(fullfile(proj,'src')));
addpath(genpath(fullfile(proj,'tools')));
rehash;

% load components and build z
S = load(fullfile(proj,'data','components.mat'),'components'); components = S.components;
z = zeros(1,numel(components)); z(1)=0.7; if numel(components)>1, z(2)=0.2; end
if numel(components)>2, z(3)=0.1; end

% baseline well and operating point
well = struct('D',0.1,'L',2000,'theta',pi/6,'roughness',1e-5);
surfaceP = 1e5; T = 320; qL = 1e-3; qG = 0.01;
opts = struct('useRefinedAnalyticJacobian',true,'refinedClampPerSeg',3e-3,'verbose',false);

% test grid
Nseg = 20;
depth_indices = [1, round(Nseg/4), round(Nseg/2), round(3*Nseg/4), Nseg]; % near-top to near-bottom
magnitudes = [2, 4, 6, 10]; % multiplicative factors on dp_seg

% prepare output CSV
folder = fileparts(outFile);
if ~exist(folder,'dir'), mkdir(folder); end
fid = fopen(outFile,'w');
fprintf(fid,'injected_segment,magnitude_factor,detected,detected_segment,severity,converged,res_norm\n');
fclose(fid);

% baseline run to get dp_seg template
res_base = NodalNewton(surfaceP, T, z, components, qL, qG, well, opts);
if ~isfield(res_base,'details') || ~isfield(res_base.details,'dp_seg')
  error('Baseline run did not return dp_seg. Ensure nodalResidualWrapper returns dp_seg.');
end
base_dp = res_base.details.dp_seg;

for s = depth_indices
  for m = magnitudes
    % copy baseline and inject synthetic leak by amplifying dp at segment s
    res = res_base;
    res.details.dp_seg = base_dp;
    res.details.dp_seg(s) = res.details.dp_seg(s) * m;
    % call detector
    report = LeakDetector(res.details, res.TubingFlow(1));
    detected = double(report.leak_flag);
    detected_segment = report.suspected_segment;
    severity = report.severity;
    converged = double(res.converged);
    res_norm = res.res_norm;
    % append row
    fid = fopen(outFile,'a');
    fprintf(fid,'%d,%d,%d,%d,%.6g,%d,%.6g\n', s, m, detected, detected_segment, severity, converged, res_norm);
    fclose(fid);
  end
end

fprintf('Validation complete. Results saved to %s\n', outFile);
end
