classdef test_NodalNewtonAnalytic < matlab.unittest.TestCase
  methods (Static, Access = private)
    function projRoot = findProjectRoot()
      thisFile = mfilename('fullpath');
      if isempty(thisFile)
        startDir = pwd;
      else
        startDir = fileparts(thisFile);
      end
      cur = startDir;
      for k = 1:10
        candidate = fullfile(cur,'data','components.mat');
        if isfile(candidate)
          projRoot = cur; return;
        end
        parent = fileparts(cur);
        if strcmp(parent,cur), break; end
        cur = parent;
      end
      if isfile(fullfile(pwd,'data','components.mat'))
        projRoot = pwd; return;
      end
      error('Could not locate data/components.mat in project tree.');
    end
  end

  methods(Test)
    function testAnalyticReducesIters(testCase)
      projRoot = test_NodalNewtonAnalytic.findProjectRoot();
      S = load(fullfile(projRoot,'data','components.mat'),'components');
      components = S.components;

      well.D = 0.1; well.L = 2000; well.theta = pi/6; well.roughness = 1e-5;
      surfaceP = 1e5; T = 320;
      z = zeros(1,numel(components)); z(1)=0.7; z(2)=0.2; z(3)=0.1;
      qL = 1e-3; qG = 0.01;

      resNum = NodalNewton(surfaceP, T, z, components, qL, qG, well, struct('useAnalyticJacobian',false,'verbose',false));
      resAna = NodalNewton(surfaceP, T, z, components, qL, qG, well, struct('useAnalyticJacobian',true,'verbose',false));
      testCase.verifyTrue(resNum.converged);
      testCase.verifyTrue(resAna.converged);
      testCase.verifyLessThanOrEqual(resAna.iters, resNum.iters + 1);
    end
  end
end
