function run_leak_detection_hysteresis(outFile, historyFile)
% run_leak_detection_hysteresis  Run cases with hysteresis on leak alarms
% Usage: run_leak_detection_hysteresis()
%        run_leak_detection_hysteresis('results/leak_reports.csv','results/leak_history.mat')

if nargin < 1 || isempty(outFile)
  outFile = fullfile(pwd,'results','leak_reports.csv');
end
if nargin < 2 || isempty(historyFile)
  historyFile = fullfile(pwd,'results','leak_history.mat');
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

% cases (same as before)
cases = {
  struct('well',struct('D',0.1,'L',2000,'theta',pi/6,'roughness',1e-5),'surfaceP',1e5,'T',320,'qL',1e-3,'qG',0.01);
  struct('well',struct('D',0.1,'L',1500,'theta',0,'roughness',1e-5),'surfaceP',1e5,'T',320,'qL',5e-4,'qG',0.005);
  struct('well',struct('D',0.08,'L',2500,'theta',pi/3,'roughness',1e-5),'surfaceP',1e5,'T',330,'qL',2e-3,'qG',0.02);
};

% prepare outputs
folder = fileparts(outFile);
if ~exist(folder,'dir'), mkdir(folder); end
if ~isfile(outFile)
  fid = fopen(outFile,'w'); fprintf(fid,'case,timestamp,leak_flag,severity,suspected_segment,BHP,res_norm,persistent_alarm\n'); fclose(fid);
end

% load history if exists
history = struct();
if isfile(historyFile)
  try load(historyFile,'history'); catch, history = struct(); end
end
if ~isfield(history,'lastSeverity'), history.lastSeverity = zeros(numel(cases),1); end
if ~isfield(history,'consecCount'), history.consecCount = zeros(numel(cases),1); end

% hysteresis parameters
persistence_required = 2; % number of consecutive runs required
severity_threshold = 0.05; % same default as LeakDetector min threshold

for i=1:numel(cases)
  c = cases{i};
  fprintf('Running case %d for leak detection with hysteresis\n', i);
  opts = struct('useRefinedAnalyticJacobian',true,'refinedClampPerSeg',3e-3,'verbose',false);
  res = NodalNewton(c.surfaceP, c.T, z, components, c.qL, c.qG, c.well, opts);

  if isfield(res,'details')
    report = LeakDetector(res.details, res.TubingFlow(1));
  else
    report = struct('leak_flag',false,'severity',0,'suspected_segment',NaN,'dp_seg',[],'p_seg',[]);
  end

  % update history counters
  if report.severity >= severity_threshold
    history.consecCount(i) = history.consecCount(i) + 1;
  else
    history.consecCount(i) = 0;
  end
  history.lastSeverity(i) = report.severity;

  % persistent alarm only if consecCount >= persistence_required
  persistent_alarm = history.consecCount(i) >= persistence_required;

  % append to CSV (include persistent_alarm flag)
  fid = fopen(outFile,'a');
  fprintf(fid,'%d,%s,%d,%.6g,%d,%.6g,%.6g,%d\n', i, char(report.timestamp), double(report.leak_flag), report.severity, report.suspected_segment, res.BHP, res.res_norm, double(persistent_alarm));
  fclose(fid);

  % save dpseg PNG as before
  try
    if ~isempty(report.dp_seg)
      Nseg = numel(report.dp_seg);
      x = linspace(0,c.well.L,Nseg+1);
      figure('Visible','off'); hold on;
      barh(x(1:end-1)+diff(x)/2, report.dp_seg, 'FaceColor',[0.6 0.6 0.9]);
      xlabel('dp per segment (Pa)'); ylabel('Depth along well (m)');
      title(sprintf('Case %d dp per segment (suspected %d)', i, report.suspected_segment));
      if ~isnan(report.suspected_segment)
        yline(x(report.suspected_segment),'r--','LineWidth',1.5);
      end
      saveas(gcf, fullfile(proj,'results',sprintf('leak_case%d_dpseg.png',i)));
      close(gcf);
    end
  catch
  end
end

% save history
save(historyFile,'history');

fprintf('Hysteresis run complete. Reports appended to %s\n', outFile);
end
