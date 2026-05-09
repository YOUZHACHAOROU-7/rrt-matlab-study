% RUN_ER7_MO_EXPERIMENT.m
% Multi-objective ER7-900 trajectory optimization launcher.
% Objective = time + jerk + energy proxy.
clear; clc; close all;
scriptFile = fullfile(pwd,'src','main_er7_mo_trajectory_experiment.m');
if ~exist(scriptFile,'file')
    error('Cannot find src/main_er7_mo_trajectory_experiment.m. Please download latest repo.');
end
run(scriptFile);
