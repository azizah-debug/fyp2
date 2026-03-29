function run_benchmarks(outFile)
% run_benchmarks  Run benchmark cases for Secant / NewtonFD / NewtonRefined
% Usage: run_benchmarks() or run_benchmarks('results/benchmarks.csv')

% Determine project root from this file's location
thisFile = mfilename('fullpath');
if isempty(thisFile)
    % fallback: assume current folder is project root
    projRoot = pwd;
else
    projRoot = fileparts(fileparts(thisFile)); % parent of tools/
end

% ensure src is on path
addpath(genpath(fullfile(projRoot,'src')));

% load components
dataFile = fullfile(projRoot,'data','components.mat');
if ~isfile(dataFile)
    error('Could not find %s. Run from project root or place components.mat in data/.', dataFile);
end
S = load(dataFile,'components'); components = S.components;

% define benchmark cases (customize as needed)
cases = {
  struct('well',struct('D',0.1,'L',2000,'theta',pi/6,'roughness',1e-5),'surfaceP',1e5,'T',320,'z',[],'qL',1e-3,'qG',0.01);
  struct('well',struct('D',0.1,'L',1500,'theta',0,'roughness',1e-5),'surfaceP',1e5,'T',320,'z',[],'qL',5e-4,'qG',0.005);
  struct('well',struct('D',0.08,'L',2500,'theta',pi/3,'roughness',1e-5),'surfaceP',1e5,'T',330,'z',[],'qL',2e-3,'qG',0.02);
};

for i=1:numel(cases)
  z = zeros(1,numel(components)); z(1)=0.7; if numel(components)>1, z(2)=0.2; end
  if numel(components)>2, z(3)=0.1; end
  cases{i}.z = z;
end

if nargin < 1 || isempty(outFile)
  outFile = fullfile(projRoot,'results','benchmarks.csv');
end

hdr = 'case,mode,converged,iters,time_s,BHP,res_norm';
folder = fileparts(outFile);
if ~exist(folder,'dir'), mkdir(folder); end
fid = fopen(outFile,'w');
fprintf(fid, '%s\n', hdr);
fclose(fid);

for ci = 1:numel(cases)
  c = cases{ci};

  % Mode A: Secant baseline (use NodalNewton FD as proxy if no secant function)
  try
    tstart = tic;
    resA = NodalNewton(c.surfaceP, c.T, c.z, components, c.qL, c.qG, c.well, struct('useAnalyticJacobian',false,'useRefinedAnalyticJacobian',false,'verbose',false));
    tA = toc(tstart);
  catch
    resA.BHP = NaN; resA.converged = false; resA.iters = NaN; resA.res_norm = NaN; tA = NaN;
  end
  writeRow(outFile, ci, 'Secant', resA, tA);

  % Mode B: Newton with FD
  try
    tstart = tic;
    resB = NodalNewton(c.surfaceP, c.T, c.z, components, c.qL, c.qG, c.well, struct('useAnalyticJacobian',false,'useRefinedAnalyticJacobian',false,'verbose',false));
    tB = toc(tstart);
  catch
    resB.BHP = NaN; resB.converged = false; resB.iters = NaN; resB.res_norm = NaN; tB = NaN;
  end
  writeRow(outFile, ci, 'NewtonFD', resB, tB);

  % Mode C: Newton with refined analytic Jacobian
  try
    tstart = tic;
    resC = NodalNewton(c.surfaceP, c.T, c.z, components, c.qL, c.qG, c.well, struct('useAnalyticJacobian',false,'useRefinedAnalyticJacobian',true,'refinedClampPerSeg',3e-3,'verbose',false));
    tC = toc(tstart);
  catch
    resC.BHP = NaN; resC.converged = false; resC.iters = NaN; resC.res_norm = NaN; tC = NaN;
  end
  writeRow(outFile, ci, 'NewtonRefined', resC, tC);
end

fprintf('Benchmarks saved to %s\n', outFile);
end

function writeRow(outFile, caseIdx, mode, res, t)
if isempty(res), res.BHP = NaN; res.converged = false; res.iters = NaN; res.res_norm = NaN; end
fid = fopen(outFile,'a');
fprintf(fid, '%d,%s,%d,%d,%.6g,%.6g,%.6g\n', caseIdx, mode, double(res.converged), res.iters, t, res.BHP, res.res_norm);
fclose(fid);
end