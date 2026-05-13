% RUN_IHHO_BENCHMARK.m
% Launcher for IHHO benchmark experiment.
% Run this file from repository root.
clear; clc; close all;

scriptFile = fullfile(pwd,'src','main_ihho_benchmark.m');
if ~exist(scriptFile,'file')
    error(['Cannot find src/main_ihho_benchmark.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end

run(scriptFile);
