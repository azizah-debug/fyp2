function compile_validation_report(validationCSV, outPDF)
% compile_validation_report  Create a multipage PDF with validation figures and summary
% Usage:
%   compile_validation_report()                     % uses defaults
%   compile_validation_report('results/validation_results.csv','results/validation_report.pdf')

if nargin < 1 || isempty(validationCSV)
  validationCSV = fullfile(pwd,'results','validation_results.csv');
end
if nargin < 2 || isempty(outPDF)
  outPDF = fullfile(pwd,'results','validation_report.pdf');
end

% Check inputs
if ~isfile(validationCSV)
  error('Validation CSV not found: %s', validationCSV);
end

% Paths to expected images
detectionRateImg = fullfile(pwd,'results','validation_detection_rate.png');
localizationHeatmapImg = fullfile(pwd,'results','validation_localization_heatmap.png');
summaryCSV = fullfile(pwd,'results','validation_summary.csv');

% Read summary if available
summaryText = {};
if isfile(summaryCSV)
  S = readtable(summaryCSV,'ReadVariableNames',false,'TextType','string');
  for i=1:height(S)
    summaryText{end+1} = sprintf('%s: %s', strtrim(char(S{i,1})), strtrim(char(S{i,2}))); %#ok<AGROW>
  end
else
  summaryText = {'No validation_summary.csv found.'};
end

% Create first page: title and summary
fig = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
axis off;
text(0.5,0.6,'Validation Report','HorizontalAlignment','center','FontSize',18,'FontWeight','bold');
nexttile;
axis off;
y = 0.9;
for k=1:numel(summaryText)
  text(0.05,y,summaryText{k},'FontSize',10,'Interpreter','none');
  y = y - 0.08;
end
drawnow;
exportgraphics(fig,outPDF,'ContentType','vector');

% Append detection rate image as second page
if isfile(detectionRateImg)
  fig2 = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
  imshow(detectionRateImg);
  title('Detection rate vs injected magnitude','Interpreter','none');
  exportgraphics(fig2,outPDF,'Append',true);
  close(fig2);
end

% Append localization heatmap as third page
if isfile(localizationHeatmapImg)
  fig3 = figure('Visible','off','Units','normalized','Position',[0 0 1 1]);
  imshow(localizationHeatmapImg);
  title('Localization error heatmap','Interpreter','none');
  exportgraphics(fig3,outPDF,'Append',true);
  close(fig3);
end

close(fig);
fprintf('Compiled validation report saved to %s\n', outPDF);
end
