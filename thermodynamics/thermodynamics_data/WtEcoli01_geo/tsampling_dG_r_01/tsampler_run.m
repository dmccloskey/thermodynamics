% initialize with Tomlab_CPLEX
load('C:\Users\dmccloskey-sbrg\Documents\MATLAB\tsampling\tsampler.mat')
initCobraToolbox();
% sample
[tsampler_out, mixedFrac] = gpSampler(iJO1366, [], [], [], [], [], true);
[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 10000, [], true);