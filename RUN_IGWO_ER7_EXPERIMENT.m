% RUN_IGWO_ER7_EXPERIMENT.m
% ER7-900 paper-level IGWO repeated experiment launcher.
clear; clc; close all;

scriptFile = fullfile(pwd,'src','main_paper_igwo_er7_experiment.m');
if ~exist(scriptFile,'file')
    error(['Cannot find src/main_paper_igwo_er7_experiment.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end

run(scriptFile);
