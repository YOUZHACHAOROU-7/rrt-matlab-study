function zipFile = finalize_experiment_results(outDir, experimentName, summaryText)
% finalize_experiment_results
% Automatically writes a summary txt file and zips the whole result folder.
%
% Usage:
%   zipFile = finalize_experiment_results(outDir,'ihho_benchmark',summaryText);
%
% Inputs:
%   outDir         result folder path, for example fullfile(pwd,'results_ihho_benchmark')
%   experimentName short experiment name used in zip filename
%   summaryText    char/string/cellstr, optional notes written to SUMMARY_README.txt
%
% Output:
%   zipFile        full path of generated zip file

    if nargin < 2 || isempty(experimentName)
        experimentName = 'experiment_results';
    end
    if nargin < 3
        summaryText = '';
    end

    if ~exist(outDir,'dir')
        error('Result folder does not exist: %s', outDir);
    end

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    summaryFile = fullfile(outDir,'SUMMARY_README.txt');

    fid = fopen(summaryFile,'w','n','UTF-8');
    if fid == -1
        warning('Cannot write SUMMARY_README.txt in %s', outDir);
    else
        fprintf(fid,'Experiment: %s\n', experimentName);
        fprintf(fid,'Generated at: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid,'Result folder: %s\n\n', outDir);
        fprintf(fid,'Files usually included:\n');
        fprintf(fid,'- *.csv: result tables for paper analysis\n');
        fprintf(fid,'- figures/*.png: convergence curves and experiment figures\n');
        fprintf(fid,'- *.mat: optional MATLAB workspace if saved by the main script\n\n');
        if iscell(summaryText)
            fprintf(fid,'Summary notes:\n');
            for i = 1:numel(summaryText)
                fprintf(fid,'- %s\n', string(summaryText{i}));
            end
        else
            fprintf(fid,'Summary notes:\n%s\n', string(summaryText));
        end
        fclose(fid);
    end

    parentDir = fileparts(outDir);
    zipFile = fullfile(parentDir, sprintf('%s_%s.zip', experimentName, timestamp));

    if exist(zipFile,'file')
        delete(zipFile);
    end

    oldDir = pwd;
    cleanupObj = onCleanup(@() cd(oldDir));
    cd(parentDir);
    [~,folderName] = fileparts(outDir);
    zip(zipFile, folderName);

    fprintf('\n========== Result Package Ready ==========\n');
    fprintf('Result folder: %s\n', outDir);
    fprintf('Zip package : %s\n', zipFile);
    fprintf('Please upload this zip package to ChatGPT for paper writing.\n');
end
