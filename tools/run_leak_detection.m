function run_leak_detection(outFile)
% run_leak_detection  Run cases and produce leak reports
% Usage: run_leak_detection() or run_leak_detection('results/leak_reports.csv')

if nargin < 1 || isempty(outFile)
  outFile = fullfile(pwd,'results','leak_reports.csv');
end

% Determine project root from this file's location
thisFile = mfilename('fullpath');
if isempty(thisFile)
    proj = pwd;
else
    proj = fileparts(fileparts(thisFile)); % parent of tools/
end
addpath(genpath(fullfile(proj,'src')));
addpath(genpath(fullfile(proj,'tools')));
rehash;


% Ensure src and tools are on the path
addpath(genpath(fullfile(proj,'src')));
addpath(genpath(fullfile(proj,'tools')));
rehash;


% load components
S = load(fullfile(proj,'data','components.mat'),'components'); components = S.components;
z = zeros(1,numel(components)); z(1)=0.7; if numel(components)>1, z(2)=0.2; end
if numel(components)>2, z(3)=0.1; end

cases = {
  struct('well',struct('D',0.1,'L',2000,'theta',pi/6,'roughness',1e-5),'surfaceP',1e5,'T',320,'qL',1e-3,'qG',0.01);
  struct('well',struct('D',0.1,'L',1500,'theta',0,'roughness',1e-5),'surfaceP',1e5,'T',320,'qL',5e-4,'qG',0.005);
  struct('well',struct('D',0.08,'L',2500,'theta',pi/3,'roughness',1e-5),'surfaceP',1e5,'T',330,'qL',2e-3,'qG',0.02);
};

% prepare output
folder = fileparts(outFile);
if ~exist(folder,'dir'), mkdir(folder); end
fid = fopen(outFile,'w');
fprintf(fid,'case,timestamp,leak_flag,severity,suspected_segment,BHP,res_norm\n');
fclose(fid);

for i=1:numel(cases)
  c = cases{i};
  fprintf('Running case %d for leak detection\n', i);
  opts = struct('useRefinedAnalyticJacobian',true,'refinedClampPerSeg',3e-3,'verbose',false);
  res = NodalNewton(c.surfaceP, c.T, z, components, c.qL, c.qG, c.well, opts);

  % call detector using details returned by NodalNewton
  if isfield(res,'details')
    report = LeakDetector(res.details, res.TubingFlow(1)); % q_in ~ liquid flow (proxy)
  else
    report = struct('leak_flag',false,'severity',0,'suspected_segment',NaN);
  end

  % append to CSV
  fid = fopen(outFile,'a');
  fprintf(fid,'%d,%s,%d,%.6g,%d,%.6g,%.6g\n', i, char(report.timestamp), double(report.leak_flag), report.severity, report.suspected_segment, res.BHP, res.res_norm);
  fclose(fid);

  % save a simple schematic PNG showing suspected segment
  try
    Nseg = numel(report.dp_seg);
    x = linspace(0,c.well.L,Nseg+1);
    figure('Visible','off'); hold on;
    barh(x(1:end-1)+diff(x)/2, report.dp_seg, 'FaceColor',[0.6 0.6 0.9]);
    xlabel('dp per segment (Pa)'); ylabel('Depth along well (m)');
    title(sprintf('Case %d dp per segment (suspected %d)', i, report.suspected_segment));
    % mark suspected segment
    if ~isnan(report.suspected_segment)
      yline(x(report.suspected_segment),'r--','LineWidth',1.5);
    end
    saveas(gcf, fullfile(proj,'results',sprintf('leak_case%d_dpseg.png',i)));
    close(gcf);
  catch
  end
end

fprintf('Leak reports saved to %s\n', outFile);
end
