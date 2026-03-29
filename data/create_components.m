% create_components.m
components = struct( ...
  'name', {'Methane','Ethane','Propane','nC4','nC5','C6plus','CO2','N2','H2O'}, ...
  'Tc',   [190.56,305.32,369.83,425.12,469.70,540.00,304.2,126.2,647.1], ... % K
  'Pc',   [4.599e6,4.884e6,4.248e6,3.800e6,3.370e6,2.800e6,7.38e6,3.39e6,22.06e6], ... % Pa
  'omega',[0.011,0.099,0.152,0.200,0.251,0.350,0.225,0.040,0.344], ...
  'Mw',   [16.04,30.07,44.10,58.12,72.15,100.0,44.01,28.01,18.02] ... % g/mol
);
save('data/components.mat','components')
disp('components.mat created')