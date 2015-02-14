% initialize with Tomlab_CPLEX
load('C:\Users\dmccloskey-sbrg\Documents\MATLAB\tsampling\tsampler.mat')
initCobraToolbox();
% sample
[tsampler_out, mixedFrac] = gpSampler(iDM2014_WTEColi_113C80_U13C20_01_OxicWtGlc, [], [], [], [], [], true);
[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 20000, [], true);