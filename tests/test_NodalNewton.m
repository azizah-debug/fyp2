classdef test_NodalNewton < matlab.unittest.TestCase
  methods(Test)
    function testNewtonConverges(testCase)
      projRoot = fileparts(fileparts(mfilename('fullpath')));
      if isempty(projRoot), projRoot = pwd; end
      S = load(fullfile(projRoot,'data','components.mat'),'components');
      components = S.components;

      well.D = 0.1; well.L = 2000; well.theta = pi/6; well.roughness = 1e-5;
      surfaceP = 1e5; T = 320;
      z = zeros(1,numel(components)); z(1)=0.7; z(2)=0.2; z(3)=0.1;
      qL = 1e-3; qG = 0.01;

      % run Newton solver
      resN = NodalNewton(surfaceP, T, z, components, qL, qG, well, struct('verbose',false));
      testCase.verifyTrue(resN.converged);
      testCase.verifyGreaterThan(resN.BHP, surfaceP);
      testCase.verifyLessThan(resN.BHP, surfaceP + 1e7);

      % compare to secant solver result (should be finite)
      resS = NodalSolverNewton(surfaceP, T, z, components, qL, qG, well);
      testCase.verifyGreaterThan(resS.BHP, surfaceP);
      testCase.verifyLessThan(resS.BHP, surfaceP + 1e7);

      % BHPs should be within a reasonable tolerance
      testCase.verifyLessThan(abs(resN.BHP - resS.BHP), 1e6);
    end
  end
end
