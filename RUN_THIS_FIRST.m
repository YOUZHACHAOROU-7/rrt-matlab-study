% RUN_THIS_FIRST.m
% Put this file in the same project root as the src folder, then click Run.
clear; clc;

if ~exist(fullfile(pwd,'src','main_ipso_353_timeopt_demo.m'),'file')
    error(['Cannot find src/main_ipso_353_timeopt_demo.m. ', ...
           'Please download the full GitHub repository or create the src folder first.']);
end

run(fullfile(pwd,'src','main_ipso_353_timeopt_demo.m'));
