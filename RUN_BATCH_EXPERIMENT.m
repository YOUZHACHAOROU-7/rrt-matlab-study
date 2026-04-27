% RUN_BATCH_EXPERIMENT.m
% Paper-level repeated experiment launcher.
% Put this file in the repository root, then click Run.
clear; clc;

batchFile = fullfile(pwd,'src','main_paper_batch_experiment.m');
if ~exist(batchFile,'file')
    error(['Cannot find src/main_paper_batch_experiment.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end

run(batchFile);
