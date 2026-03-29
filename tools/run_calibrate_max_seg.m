function run_calibrate_max_seg()
% run_calibrate_max_seg  calibrate per-segment clamp for refined Jacobian
% Saves results to results/calibration_max_seg.mat

% simple calibration set (tune these cases to match your project)
projRoot = fileparts(fileparts(mfilename('fullpath')));
dataFile = fullfile(projRoot,'data','components.mat');
if ~isfile(dataFile)
    error('Could not find data/components.mat. Run from project root or place components.mat in data/.');
end
S = load(dataFile,'components'); components = S.components;

% define a few representative wells
cases = {
    struct('well',struct('D',0.1,'L',2000,'theta',pi/6,'roughness',1e-5),'surfaceP',1e5,'T',320,'z',[],'qL',1e-3,'qG',0.01);
    struct('well',struct('D',0.1,'L',1500,'theta',0,'roughness',1e-5),'surfaceP',1e5,'T',320,'z',[],'qL',5e-4,'qG',0.005);
    struct('well',struct('D',0.08,'L',2500,'theta',pi/3,'roughness',1e-5),'surfaceP',1e5,'T',330,'z',[],'qL',2e-3,'qG',0.02);
};

% default composition vector (adjust to your components)
for i=1:numel(cases)
    z = zeros(1,numel(components)); z(1)=0.7; z(2)=0.2; if numel(components)>2, z(3)=0.1; end
    cases{i}.z = z;
end

% clamp candidates to try
clamps = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1];

results = struct();
for ci = 1:numel(clamps)
    perSegClamp = clamps(ci);
    iters = zeros(1,numel(cases));
    times = zeros(1,numel(cases));
    for k=1:numel(cases)
        c = cases{k};
        opts = struct('useRefinedAnalyticJacobian',true,'useAnalyticJacobian',false,'verbose',false,'refinedClampPerSeg',perSegClamp);
        tic;
        res = NodalNewton(c.surfaceP, c.T, c.z, components, c.qL, c.qG, c.well, opts);
        times(k) = toc;
        iters(k) = res.iters;
    end
    results(ci).perSegClamp = perSegClamp;
    results(ci).iters = iters;
    results(ci).times = times;
    results(ci).meanIter = mean(iters);
    results(ci).medianTime = median(times);
    fprintf('Clamp %g -> mean iters %.2f, median time %.3fs\n', perSegClamp, results(ci).meanIter, results(ci).medianTime);
end

outFile = fullfile(projRoot,'results','calibration_max_seg.mat');
if ~exist(fileparts(outFile),'dir'), mkdir(fileparts(outFile)); end
save(outFile,'results','clamps');
fprintf('Calibration saved to %s\n', outFile);
end
