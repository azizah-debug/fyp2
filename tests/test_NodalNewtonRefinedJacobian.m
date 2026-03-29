classdef test_NodalNewtonRefinedJacobian < matlab.unittest.TestCase
  methods(Test)
    function testRefinedHelps(testCase)
      projRoot = fileparts(fileparts(mfilename('fullpath')));
      if isempty(projRoot), projRoot = pwd; end
      S = load(fullfile(projRoot,'data','components.mat'),'components');
      components = S.components;

      well.D = 0.1; well.L = 2000; well.theta = pi/6; well.roughness = 1e-5;
      surfaceP = 1e5; T = 320;
      z = zeros(1,numel(components)); z(1)=0.7; z(2)=0.2; z(3)=0.1;
      qL = 1e-3; qG = 0.01;

      resNum = NodalNewton(surfaceP, T, z, components, qL, qG, well, struct('useAnalyticJacobian',false,'useRefinedAnalyticJacobian',false,'verbose',false));
      resRef = NodalNewton(surfaceP, T, z, components, qL, qG, well, struct('useAnalyticJacobian',false,'useRefinedAnalyticJacobian',true,'verbose',false));
      testCase.verifyTrue(resNum.converged);
      testCase.verifyTrue(resRef.converged);
      testCase.verifyLessThanOrEqual(resRef.iters, resNum.iters + 1);
    end
  end
end
