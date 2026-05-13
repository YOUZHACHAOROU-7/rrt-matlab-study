% RUN_ALGORITHM_SCREENING.m
% Candidate algorithm screening launcher for the ER7-900 trajectory planning paper.
% It compares IGWO/HGWO/TFHO/IPSO and several baseline algorithms on benchmark functions.
clear; clc; close all;

scriptFile = fullfile(pwd,'src','main_algorithm_screening.m');
if ~exist(scriptFile,'file')
    error(['Cannot find src/main_algorithm_screening.m. ', ...
           'Please download the latest GitHub repository main branch.']);
end

run(scriptFile);
