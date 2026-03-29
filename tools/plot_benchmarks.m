function plot_benchmarks(csvFile)
% plot_benchmarks  Plot benchmark iterations and runtime
% Usage: plot_benchmarks() or plot_benchmarks('results/benchmarks.csv')

% Determine project root from this file's location
thisFile = mfilename('fullpath');
if isempty(thisFile)
    projRoot = pwd;
else
    projRoot = fileparts(fileparts(thisFile)); % parent of tools/
end

if nargin < 1 || isempty(csvFile)
    csvFile = fullfile(projRoot,'results','benchmarks.csv');
end

if ~isfile(csvFile)
    error('Could not find %s. Make sure run_benchmarks saved the file.', csvFile);
end

% Read table, preserve original headers if present
T = readtable(csvFile,'PreserveVariableNames',true);

% Normalize variable names if they were modified
if any(strcmpi(T.Properties.VariableNames,'xCase'))
    % older output used xCase; rename for convenience
    T.Properties.VariableNames = {'case','mode','converged','iters','time_s','BHP','res_norm'};
end

% Convert mode to categorical for plotting order
modes = unique(T.mode,'stable');

% Plot iterations
figure('Name','Iterations by Mode','NumberTitle','off');
hold on;
for i=1:numel(modes)
    rows = strcmp(T.mode,modes{i});
    plot(find(rows), T.iters(rows), '-o','DisplayName',modes{i});
end
xlabel('Run index');
ylabel('Iterations');
legend('Location','best');
grid on;
saveas(gcf, fullfile(projRoot,'results','bench_iters.png'));

% Plot time
figure('Name','Time by Mode','NumberTitle','off');
hold on;
for i=1:numel(modes)
    rows = strcmp(T.mode,modes{i});
    plot(find(rows), T.time_s(rows), '-o','DisplayName',modes{i});
end
xlabel('Run index');
ylabel('Time (s)');
legend('Location','best');
grid on;
saveas(gcf, fullfile(projRoot,'results','bench_time.png'));

fprintf('Plots saved to %s and %s\n', fullfile(projRoot,'results','bench_iters.png'), fullfile(projRoot,'results','bench_time.png'));
end
