classdef test_PR_unitTests < matlab.unittest.TestCase
  methods (Static, Access = private)
    function projRoot = findProjectRoot()
      % Start from this test file location and walk up until data/components.mat is found
      thisFile = mfilename('fullpath');
      if isempty(thisFile)
        % fallback to current folder
        startDir = pwd;
      else
        startDir = fileparts(thisFile);
      end
      cur = startDir;
      maxDepth = 10;
      found = false;
      for k = 1:maxDepth
        candidate = fullfile(cur,'data','components.mat');
        if isfile(candidate)
          projRoot = cur;
          found = true;
          return;
        end
        parent = fileparts(cur);
        if strcmp(parent,cur)
          break;
        end
        cur = parent;
      end
      if ~found
        % final fallback: try current working folder
        if isfile(fullfile(pwd,'data','components.mat'))
          projRoot = pwd;
          return;
        end
        error('Could not locate data/components.mat in project tree. Ensure components.mat exists under a data folder.');
      end
    end
  end

  methods(Test)
    function testMethaneZ(testCase)
      projRoot = test_PR_unitTests.findProjectRoot();
      compPath = fullfile(projRoot,'data','components.mat');
      S = load(compPath,'components');
      components = S.components;

      x = zeros(1,numel(components)); x(1)=1;
      Z = PR_eos(1e7,300,components,x,'vapor');
      testCase.verifyGreaterThan(Z,0.2);
      testCase.verifyLessThan(Z,1.5);
    end

    function testPRMixPropsDensity(testCase)
      projRoot = test_PR_unitTests.findProjectRoot();
      compPath = fullfile(projRoot,'data','components.mat');
      S = load(compPath,'components');
      components = S.components;

      z = zeros(1,numel(components)); z(1)=1;
      props = PR_mixProps(1e7,300,components,z,'vapor');
      testCase.verifyGreaterThan(props.rho, 0.1);
      testCase.verifyLessThan(props.rho, 1000);
    end

    function testIsothermalFlashMassBalance(testCase)
      projRoot = test_PR_unitTests.findProjectRoot();
      compPath = fullfile(projRoot,'data','components.mat');
      S = load(compPath,'components');
      components = S.components;

      z = zeros(1,numel(components)); z(1)=0.7; z(2)=0.2; z(3)=0.1;
      [beta, x, y] = isothermalFlash(5e6, 320, components, z);
      z_rec = (1-beta)*x + beta*y;
      testCase.verifyLessThan(max(abs(z_rec - z)), 1e-6);
    end
  end
end