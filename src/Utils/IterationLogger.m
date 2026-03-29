classdef IterationLogger < handle
  properties
    fid = -1;
    filename = '';
    startTime = 0;
  end

  methods
    function obj = IterationLogger(outFile)
      if nargin < 1 || isempty(outFile)
        outFile = fullfile(pwd,'results','iteration_log.csv');
      end
      obj.filename = outFile;
      folder = fileparts(outFile);
      if ~exist(folder,'dir'), mkdir(folder); end
      obj.fid = fopen(outFile,'w');
      if obj.fid < 0
        error('Could not open log file %s for writing', outFile);
      end
      fprintf(obj.fid, 'timestamp,iter,BHP,abs_res,jac_source,alpha,dp,elapsed_s\n');
      obj.startTime = tic;
    end

    function log(obj, iter, BHP, res, jacSource, alpha, dp)
      if obj.fid < 0, return; end
      ts = datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
      elapsed = toc(obj.startTime);
      if isempty(alpha), alpha = NaN; end
      if isempty(dp), dp = NaN; end
      fprintf(obj.fid, '%s,%d,%.6g,%.6g,%s,%.6g,%.6g,%.6g\n', ...
        char(ts), iter, BHP, abs(res), jacSource, alpha, dp, elapsed);
      fflush(obj.fid);
    end

    function close(obj)
      if obj.fid > 0
        fclose(obj.fid);
        obj.fid = -1;
      end
    end
  end
end
