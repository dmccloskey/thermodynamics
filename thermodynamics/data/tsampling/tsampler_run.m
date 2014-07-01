% initialize with Tomlab_CPLEX
load('C:\Users\dmccloskey-sbrg\Documents\MATLAB	sampling	sampler_run.mat')
initCobraToolbox();
load('iJO1366');
% sample
[tsampler_out, mixedFrac] = gpSampler(tsampler, [], [], [], [], [], true);
[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 10000, [], true);