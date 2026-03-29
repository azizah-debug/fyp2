classdef test_BeggsBrill < matlab.unittest.TestCase
  methods(Test)
    function testHorizontalCase(testCase)
      out = BeggsBrill_full(1e-3, 0.01, 0.1, 0, 800, 50, 1e-3, 1e-5, 0.02);
      testCase.verifyGreaterThan(out.dp_dz, 0);
      testCase.verifyGreaterThan(out.HL, 0);
    end
    function testInclinedCase(testCase)
      out = BeggsBrill_full(5e-4, 0.005, 0.1, pi/6, 850, 60, 1e-3, 1e-5, 0.02);
      testCase.verifyLessThan(out.dp_dz, 1e5);
    end
  end
end
