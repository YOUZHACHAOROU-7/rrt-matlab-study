% RUN_BENCHMARK_IGWO.m
% Benchmark launcher for paper-level IGWO comparison.
% Run this file from the repository root.
clear; clc; close all;

scriptFile = fullfile(pwd,'src','main_benchmark_igwo_compare.m');
if ~exist(scriptFile,'file')
    error(['Cannot find src/main_benchmark_igwo_compare.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end

run(scriptFile);
