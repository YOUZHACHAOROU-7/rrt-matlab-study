% RUN_ER7_MULTI_OBJECTIVE_EXPERIMENT.m
% ER7-900 multi-objective trajectory optimization launcher.
% Objective: time + jerk + energy proxy.
clear; clc; close all;
scriptFile = fullfile(pwd,'src','main_er7_multiobjective_experiment.m');
if ~exist(scriptFile,'file')
    error(['Cannot find src/main_er7_multiobjective_experiment.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end
run(scriptFile);
