%% main_ihho_tuning.m
% Quick tuning/check script for IHHO stability on difficult functions.
% Runs only F4, F8 and F14 with fewer repeated runs.

clear; clc; close all; rng(20260513);

rootDir = pwd;
if endsWith(rootDir,[filesep 'src'])
    rootDir = fileparts(rootDir);
end
addpath(fullfile(rootDir,'src'));

outDir = fullfile(rootDir,'results_ihho_tuning');
if ~exist(outDir,'dir'), mkdir(outDir); end
figDir = fullfile(outDir,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

N = 30;
MaxIt = 300;
RunNo = 10;
Dim = 30;
algNames = {'HHO','IHHO'};
funcNames = {'F4_Rosenbrock','F8_Levy','F14_Schwefel226'};

Best = zeros(numel(funcNames),numel(algNames));
Mean = Best; Std = Best;

fprintf('IHHO quick tuning starts: %d functions, %d runs.\n',numel(funcNames),RunNo);

for f = 1:numel(funcNames)
    [lb,ub,dim,fobj] = getBenchmark(funcNames{f},Dim);
    bestRun = zeros(RunNo,numel(algNames));
    curvesRun = zeros(MaxIt,numel(algNames),RunNo);
    fprintf('\n==== %s ====\n',funcNames{f});
    for r = 1:RunNo
        rng(900000 + f*1000 + r);
        [bestRun(r,1),~,curvesRun(:,1,r)] = algHHO(N,MaxIt,lb,ub,dim,fobj);
        [bestRun(r,2),~,curvesRun(:,2,r)] = IHHO(N,MaxIt,lb,ub,dim,fobj);
        fprintf('Run %02d/%02d | HHO %.3e | IHHO %.3e\n',r,RunNo,bestRun(r,1),bestRun(r,2));
    end
    Best(f,:) = min(bestRun,[],1);
    Mean(f,:) = mean(bestRun,1);
    Std(f,:) = std(bestRun,0,1);
    T = table((1:RunNo)',bestRun(:,1),bestRun(:,2),'VariableNames',{'RunID','HHO','IHHO'});
    writetable(T,fullfile(outDir,[funcNames{f} '_runs.csv']));

    avgCurve = mean(curvesRun,3);
    fig = figure('Color','w','Visible','off');
    semilogy(abs(avgCurve(:,1))+eps,'LineWidth',1.5); hold on;
    semilogy(abs(avgCurve(:,2))+eps,'LineWidth',1.8);
    grid on; xlabel('Iteration'); ylabel('Mean best fitness');
    title(['HHO vs IHHO on ' strrep(funcNames{f},'_','\_')]);
    legend('HHO','IHHO','Location','northeast');
    saveas(fig,fullfile(figDir,[funcNames{f} '_tuning_convergence.png']));
    close(fig);
end

Summary = table(funcNames',Best(:,1),Mean(:,1),Std(:,1),Best(:,2),Mean(:,2),Std(:,2), ...
    'VariableNames',{'Function','HHO_Best','HHO_Mean','HHO_Std','IHHO_Best','IHHO_Mean','IHHO_Std'});
writetable(Summary,fullfile(outDir,'tuning_summary.csv'));
disp('========== IHHO Quick Tuning Summary ==========');
disp(Summary);

summaryText = {
    'IHHO quick tuning experiment finished.';
    'This script checks IHHO stability on F4, F8 and F14.';
    'Upload the generated zip file to ChatGPT if further tuning is needed.'
};
if exist('finalize_experiment_results','file') == 2
    finalize_experiment_results(outDir,'ihho_tuning',summaryText);
end

%% local functions copied for standalone tuning
function [lb,ub,dim,fobj] = getBenchmark(name,dim)
switch name
    case 'F4_Rosenbrock'
        lb=-30*ones(1,dim); ub=30*ones(1,dim); fobj=@(x)sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
    case 'F8_Levy'
        lb=-10*ones(1,dim); ub=10*ones(1,dim); fobj=@(x)levyFun(x);
    case 'F14_Schwefel226'
        lb=-500*ones(1,dim); ub=500*ones(1,dim); fobj=@(x)418.9829*numel(x)-sum(x.*sin(sqrt(abs(x))));
    otherwise
        error('Unknown function %s',name);
end
end
function y = levyFun(x)
w = 1 + (x-1)/4;
y = sin(pi*w(1))^2 + sum((w(1:end-1)-1).^2 .* (1+10*sin(pi*w(1:end-1)+1).^2)) + (w(end)-1)^2 * (1+sin(2*pi*w(end))^2);
end
function X = initRand(N,dim,lb,ub), X = rand(N,dim).*(ub-lb)+lb; end
function X = bound(X,lb,ub), X = min(max(X,lb),ub); end
function step = levyStep(dim,beta)
sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = randn(1,dim)*sigma; v = randn(1,dim); step = u./(abs(v).^(1/beta)+eps); step = max(min(step,10),-10);
end
function [bestF,bestX,curve] = algHHO(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); bestF = inf; bestX = zeros(1,dim); curve = zeros(MaxIt,1);
for it = 1:MaxIt
    for i = 1:N
        X(i,:) = bound(X(i,:),lb,ub); f = fobj(X(i,:));
        if f < bestF, bestF = f; bestX = X(i,:); end
    end
    E1 = 2*(1-it/MaxIt);
    for i = 1:N
        E0 = 2*rand-1; E = E1*E0; q = rand; r = rand; J = 2*(1-rand);
        if abs(E) >= 1
            randX = X(randi(N),:);
            if q < 0.5
                X(i,:) = randX - rand(1,dim).*abs(randX-2*rand(1,dim).*X(i,:));
            else
                X(i,:) = (bestX-mean(X,1)) - rand(1,dim).*((ub-lb).*rand(1,dim)+lb);
            end
        else
            if r >= 0.5 && abs(E) < 0.5
                X(i,:) = bestX - E.*abs(bestX-X(i,:));
            elseif r >= 0.5 && abs(E) >= 0.5
                X(i,:) = (bestX-X(i,:)) - E.*abs(J*bestX-X(i,:));
            else
                Y = bestX - E.*abs(J*bestX-X(i,:));
                Z = Y + randn(1,dim).*levyStep(dim,1.5);
                if fobj(bound(Y,lb,ub)) < fobj(X(i,:)), X(i,:) = Y;
                elseif fobj(bound(Z,lb,ub)) < fobj(X(i,:)), X(i,:) = Z; end
            end
        end
    end
    curve(it) = bestF;
end
end
