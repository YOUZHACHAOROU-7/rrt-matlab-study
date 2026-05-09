% RUN_ER7_MO_TIMEPRIORITY.m
% Time-priority multi-objective ER7-900 trajectory optimization.
% This runner creates a temporary copy of src/main_er7_mo_trajectory_experiment.m
% and changes only the objective weights and output folder.
%
% Objective weights:
%   Time  = 0.80
%   Jerk  = 0.10
%   Energy= 0.10
%
% Use this result as the main ER7 experiment if it shortens time while keeping jerk/energy acceptable.

clear; clc; close all;
rootDir = pwd;
if endsWith(rootDir,[filesep 'src'])
    rootDir = fileparts(rootDir);
end
srcFile = fullfile(rootDir,'src','main_er7_mo_trajectory_experiment.m');
if ~exist(srcFile,'file')
    error('Cannot find src/main_er7_mo_trajectory_experiment.m. Please download latest repo first.');
end
code = fileread(srcFile);
code = strrep(code, "results_er7_multiobjective", "results_er7_mo_timepriority");
code = strrep(code, "opt.wTime = 0.50; opt.wJerk = 0.25; opt.wEnergy = 0.25;", ...
                   "opt.wTime = 0.80; opt.wJerk = 0.10; opt.wEnergy = 0.10;");
code = strrep(code, "ER7 multi-objective experiment starts", "ER7 time-priority multi-objective experiment starts");

tmpFile = fullfile(rootDir,'src','__tmp_er7_mo_timepriority_experiment.m');
fid = fopen(tmpFile,'w');
fwrite(fid,code);
fclose(fid);
run(tmpFile);
