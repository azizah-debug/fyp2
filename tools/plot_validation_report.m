function plot_validation_report(validationCSV)
% plot_validation_report  Read validation CSV, compute metrics, and save plots
% Usage: plot_validation_report('results/validation_results.csv')

if nargin < 1 || isempty(validationCSV)
  validationCSV = fullfile(pwd,'results','validation_results.csv');
end

if ~isfile(validationCSV)
  error('Validation CSV not found: %s', validationCSV);
end

T = readtable(validationCSV,'PreserveVariableNames',true);

% Ensure expected columns exist
req = {'injected_segment','magnitude_factor','detected','detected_segment','severity','converged','res_norm'};
for k=1:numel(req)
  if ~ismember(req{k}, T.Properties.VariableNames)
    error('Missing column %s in %s', req{k}, validationCSV);
  end
end

% Convert to numeric arrays
injected = T.injected_segment;
mags = T.magnitude_factor;
detected = T.detected;
detected_seg = T.detected_segment;
severity = T.severity;
converged = T.converged;

% Unique magnitudes and segments
unique_mags = unique(mags);
unique_segs = unique(injected);

% Detection rate vs magnitude
det_rate = zeros(size(unique_mags));
for i=1:numel(unique_mags)
  idx = mags == unique_mags(i);
  det_rate(i) = sum(detected(idx)) / sum(idx);
end

% Plot detection rate
figure('Visible','off'); hold on;
plot(unique_mags, det_rate, '-o','LineWidth',1.5);
xlabel('Magnitude factor'); ylabel('Detection rate');
title('Detection rate vs injected magnitude');
grid on;
saveas(gcf, fullfile(pwd,'results','validation_detection_rate.png'));
close(gcf);

% Localization error matrix
err_mat = nan(numel(unique_segs), numel(unique_mags));
det_rate_mat = nan(numel(unique_segs), numel(unique_mags));
for i=1:numel(unique_segs)
  for j=1:numel(unique_mags)
    idx = injected==unique_segs(i) & mags==unique_mags(j);
    if any(idx)
      detected_idx = detected(idx) == 1;
      if any(detected_idx)
        err = abs(detected_seg(idx(detected_idx)) - injected(idx(detected_idx)));
        err_mat(i,j) = mean(err);
      else
        err_mat(i,j) = NaN;
      end
      det_rate_mat(i,j) = sum(detected(idx)) / sum(idx);
    end
  end
end

% Heatmap of localization error
figure('Visible','off');
imagesc(unique_mags, unique_segs, err_mat);
colorbar;
xlabel('Magnitude factor'); ylabel('Injected segment index');
title('Mean localization error');
set(gca,'YDir','normal');
saveas(gcf, fullfile(pwd,'results','validation_localization_heatmap.png'));
close(gcf);

% Summary metrics
total_tests = height(T);
overall_detection_rate = sum(detected) / total_tests;

% Save a small summary CSV robustly: write to temp then move
resultsDir = fullfile(pwd,'results');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

tmp = [tempname '.csv'];
fid = fopen(tmp,'w');
if fid == -1
  error('plot_validation_report:FileOpenError','Cannot open temp summary file for writing: %s', tmp);
end
fprintf(fid,'metric,value\n');
fprintf(fid,'total_tests,%d\n', total_tests);
fprintf(fid,'overall_detection_rate,%.6g\n', overall_detection_rate);
fclose(fid);

summaryFile = fullfile(resultsDir,'validation_summary.csv');
try
  % If a file with the same name exists, attempt to move it aside first
  if isfile(summaryFile)
    backup = fullfile(resultsDir, ['validation_summary_backup_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
    movefile(summaryFile, backup);
  end
  movefile(tmp, summaryFile);
catch ME
  % If move fails, leave temp file and report the path
  warning('plot_validation_report:MoveFailed','Could not move temp summary to results: %s\nTemp file left at: %s', ME.message, tmp);
  return;
end

fprintf('Validation plots saved to results/validation_detection_rate.png and results/validation_localization_heatmap.png\n');
fprintf('Summary saved to %s\n', summaryFile);
end
