%% main_ihho_benchmark.m
% Benchmark experiment for paper-level IHHO verification.
% Algorithms: PSO, GWO, WOA, DE, SSA, HHO and IHHO.
% Outputs: CSV result tables, PNG convergence curves, MAT workspace and ZIP package.

clear; clc; close all; rng(20260513);

rootDir = pwd;
if endsWith(rootDir,[filesep 'src'])
    rootDir = fileparts(rootDir);
end
addpath(fullfile(rootDir,'src'));

outDir = fullfile(rootDir,'results_ihho_benchmark');
if ~exist(outDir,'dir'), mkdir(outDir); end
figDir = fullfile(outDir,'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

if exist('IHHO','file') ~= 2
    error('Cannot find src/IHHO.m. Please download the latest repository main branch.');
end

%% Settings
N = 30;          % population size
MaxIt = 500;     % max iterations
RunNo = 30;      % independent runs; set to 3 for quick debugging
Dim = 30;        % benchmark dimension

algNames = {'PSO','GWO','WOA','DE','SSA','HHO','IHHO'};
funcNames = {'F1_Sphere','F2_Schwefel222','F3_Schwefel12','F4_Rosenbrock',...
             'F5_Rastrigin','F6_Ackley','F7_Griewank','F8_Levy',...
             'F9_Alpine','F10_QuarticNoise','F11_Zakharov','F12_Salomon',...
             'F13_Step','F14_Schwefel226'};

nAlg = numel(algNames);
nFunc = numel(funcNames);
Best = zeros(nFunc,nAlg);
Mean = zeros(nFunc,nAlg);
Std = zeros(nFunc,nAlg);
Median = zeros(nFunc,nAlg);
Rank = zeros(nFunc,nAlg);
FirstFlag = zeros(nFunc,nAlg);
Top2Flag = zeros(nFunc,nAlg);
AllBest = zeros(RunNo,nAlg,nFunc);
AllCurves = zeros(MaxIt,nAlg,nFunc);

fprintf('IHHO benchmark starts: %d functions, %d algorithms, %d runs.\n',nFunc,nAlg,RunNo);

for f = 1:nFunc
    [lb,ub,dim,fobj] = getBenchmark(funcNames{f},Dim);
    curvesRun = zeros(MaxIt,nAlg,RunNo);
    bestRun = zeros(RunNo,nAlg);
    fprintf('\n==== %s ====\n',funcNames{f});

    for r = 1:RunNo
        rng(20260513 + f*1000 + r);
        for a = 1:nAlg
            [score,~,curve] = runAlgorithm(algNames{a},N,MaxIt,lb,ub,dim,fobj);
            bestRun(r,a) = score;
            curvesRun(:,a,r) = curve(:);
        end
        fprintf('Run %02d/%02d | ',r,RunNo);
        for a = 1:nAlg
            fprintf('%s %.3e ',algNames{a},bestRun(r,a));
        end
        fprintf('\n');
    end

    AllBest(:,:,f) = bestRun;
    avgCurve = mean(curvesRun,3);
    AllCurves(:,:,f) = avgCurve;
    Best(f,:) = min(bestRun,[],1);
    Mean(f,:) = mean(bestRun,1);
    Std(f,:) = std(bestRun,0,1);
    Median(f,:) = median(bestRun,1);

    [~,ord] = sort(Mean(f,:),'ascend');
    for k = 1:nAlg
        Rank(f,ord(k)) = k;
    end
    FirstFlag(f,ord(1)) = 1;
    Top2Flag(f,ord(1:2)) = 1;

    T = table((1:RunNo)','VariableNames',{'RunID'});
    for a = 1:nAlg
        T.(algNames{a}) = bestRun(:,a);
    end
    writetable(T,fullfile(outDir,[funcNames{f} '_runs.csv']));

    fig = figure('Color','w','Name',funcNames{f},'Visible','off'); hold on;
    for a = 1:nAlg
        y = avgCurve(:,a);
        valid = isfinite(y);
        if ~any(valid)
            y = ones(size(y));
        else
            y(~valid) = max(y(valid));
        end
        if all(y <= 0)
            plot(y,'LineWidth',1.4);
        else
            semilogy(abs(y)+eps,'LineWidth',1.4);
        end
    end
    grid on; xlabel('Iteration'); ylabel('Mean best fitness');
    title(['Average convergence curves on ' strrep(funcNames{f},'_','\_')]);
    legend(algNames,'Location','northeastoutside');
    saveas(fig,fullfile(figDir,[funcNames{f} '_convergence.png']));
    close(fig);
end

%% Summary tables
Summary = table(funcNames','VariableNames',{'Function'});
for a = 1:nAlg
    Summary.([algNames{a} '_Best']) = Best(:,a);
    Summary.([algNames{a} '_Mean']) = Mean(:,a);
    Summary.([algNames{a} '_Std']) = Std(:,a);
    Summary.([algNames{a} '_Median']) = Median(:,a);
    Summary.([algNames{a} '_Rank']) = Rank(:,a);
end
writetable(Summary,fullfile(outDir,'ihho_benchmark_summary.csv'));

AvgRank = mean(Rank,1)';
FirstCount = sum(FirstFlag,1)';
Top2Count = sum(Top2Flag,1)';
MeanOfMean = mean(Mean,1)';
RankTable = table(algNames',AvgRank,FirstCount,Top2Count,MeanOfMean, ...
    'VariableNames',{'Algorithm','AverageRank','FirstCount','Top2Count','MeanOfMeanFitness'});
RankTable = sortrows(RankTable,'AverageRank','ascend');
writetable(RankTable,fullfile(outDir,'ihho_benchmark_rank_table.csv'));

save(fullfile(outDir,'ihho_benchmark_workspace.mat'), ...
    'algNames','funcNames','Best','Mean','Std','Median','Rank','FirstFlag','Top2Flag','AllBest','AllCurves','RankTable');

disp('========== IHHO Benchmark Rank Table ==========');
disp(RankTable);
fprintf('\nAll IHHO benchmark tables and figures saved to: %s\n',outDir);

summaryText = {
    'This experiment compares PSO, GWO, WOA, DE, SSA, HHO and IHHO on 14 benchmark functions.';
    'Each function is tested with repeated independent runs under identical population size and iteration settings.';
    'CSV tables, PNG convergence curves and MAT workspace have been saved automatically.';
    'Use ihho_benchmark_rank_table.csv and ihho_benchmark_summary.csv for paper result analysis.'
};

if exist('finalize_experiment_results','file') == 2
    finalize_experiment_results(outDir,'ihho_benchmark',summaryText);
else
    warning('finalize_experiment_results.m was not found. ZIP package was not generated.');
end

%% ====================== local functions ======================
function [score,pos,curve] = runAlgorithm(name,N,MaxIt,lb,ub,dim,fobj)
switch upper(name)
    case 'PSO'
        [score,pos,curve] = algPSO(N,MaxIt,lb,ub,dim,fobj);
    case 'GWO'
        [score,pos,curve] = algGWO(N,MaxIt,lb,ub,dim,fobj);
    case 'WOA'
        [score,pos,curve] = algWOA(N,MaxIt,lb,ub,dim,fobj);
    case 'DE'
        [score,pos,curve] = algDE(N,MaxIt,lb,ub,dim,fobj);
    case 'SSA'
        [score,pos,curve] = algSSA(N,MaxIt,lb,ub,dim,fobj);
    case 'HHO'
        [score,pos,curve] = algHHO(N,MaxIt,lb,ub,dim,fobj);
    case 'IHHO'
        [score,pos,curve] = IHHO(N,MaxIt,lb,ub,dim,fobj);
    otherwise
        error('Unknown algorithm: %s',name);
end
end

function [lb,ub,dim,fobj] = getBenchmark(name,dim)
switch name
    case 'F1_Sphere'
        lb=-100*ones(1,dim); ub=100*ones(1,dim); fobj=@(x)sum(x.^2);
    case 'F2_Schwefel222'
        lb=-10*ones(1,dim); ub=10*ones(1,dim); fobj=@(x)sum(abs(x))+prod(abs(x));
    case 'F3_Schwefel12'
        lb=-100*ones(1,dim); ub=100*ones(1,dim); fobj=@(x)sum(cumsum(x).^2);
    case 'F4_Rosenbrock'
        lb=-30*ones(1,dim); ub=30*ones(1,dim); fobj=@(x)sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
    case 'F5_Rastrigin'
        lb=-5.12*ones(1,dim); ub=5.12*ones(1,dim); fobj=@(x)sum(x.^2-10*cos(2*pi*x)+10);
    case 'F6_Ackley'
        lb=-32*ones(1,dim); ub=32*ones(1,dim); fobj=@(x)-20*exp(-0.2*sqrt(sum(x.^2)/numel(x)))-exp(sum(cos(2*pi*x))/numel(x))+20+exp(1);
    case 'F7_Griewank'
        lb=-600*ones(1,dim); ub=600*ones(1,dim); fobj=@(x)sum(x.^2)/4000-prod(cos(x./sqrt(1:numel(x))))+1;
    case 'F8_Levy'
        lb=-10*ones(1,dim); ub=10*ones(1,dim); fobj=@(x)levyFun(x);
    case 'F9_Alpine'
        lb=-10*ones(1,dim); ub=10*ones(1,dim); fobj=@(x)sum(abs(x.*sin(x)+0.1*x));
    case 'F10_QuarticNoise'
        lb=-1.28*ones(1,dim); ub=1.28*ones(1,dim); fobj=@(x)sum((1:numel(x)).*(x.^4))+rand;
    case 'F11_Zakharov'
        lb=-5*ones(1,dim); ub=10*ones(1,dim); fobj=@(x)sum(x.^2)+(sum(0.5*(1:numel(x)).*x))^2+(sum(0.5*(1:numel(x)).*x))^4;
    case 'F12_Salomon'
        lb=-100*ones(1,dim); ub=100*ones(1,dim); fobj=@(x)1-cos(2*pi*sqrt(sum(x.^2)))+0.1*sqrt(sum(x.^2));
    case 'F13_Step'
        lb=-100*ones(1,dim); ub=100*ones(1,dim); fobj=@(x)sum(floor(x+0.5).^2);
    case 'F14_Schwefel226'
        lb=-500*ones(1,dim); ub=500*ones(1,dim); fobj=@(x)418.9829*numel(x)-sum(x.*sin(sqrt(abs(x))));
    otherwise
        error('Unknown benchmark function: %s',name);
end
end

function y = levyFun(x)
w = 1 + (x-1)/4;
y = sin(pi*w(1))^2 + sum((w(1:end-1)-1).^2 .* (1+10*sin(pi*w(1:end-1)+1).^2)) + ...
    (w(end)-1)^2 * (1+sin(2*pi*w(end))^2);
end

function X = initRand(N,dim,lb,ub)
X = rand(N,dim).*(ub-lb)+lb;
end

function X = bound(X,lb,ub)
X = min(max(X,lb),ub);
end

function step = levyStep(dim,beta)
sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = randn(1,dim)*sigma;
v = randn(1,dim);
step = u./(abs(v).^(1/beta)+eps);
step = max(min(step,10),-10);
end

function [bestF,bestX,curve] = algPSO(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); V = zeros(N,dim); P = X; Pf = inf(N,1); curve = zeros(MaxIt,1);
for i = 1:N, Pf(i) = fobj(X(i,:)); end
[bestF,id] = min(Pf); bestX = P(id,:);
for it = 1:MaxIt
    w = 0.9 - 0.5*it/MaxIt; c1 = 2; c2 = 2;
    for i = 1:N
        V(i,:) = w*V(i,:) + c1*rand(1,dim).*(P(i,:)-X(i,:)) + c2*rand(1,dim).*(bestX-X(i,:));
        X(i,:) = bound(X(i,:)+V(i,:),lb,ub);
        f = fobj(X(i,:));
        if f < Pf(i), Pf(i) = f; P(i,:) = X(i,:); end
        if f < bestF, bestF = f; bestX = X(i,:); end
    end
    curve(it) = bestF;
end
end

function [Alpha_score,Alpha_pos,curve] = algGWO(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); curve = zeros(MaxIt,1);
Alpha_score = inf; Beta_score = inf; Delta_score = inf;
Alpha_pos = zeros(1,dim); Beta_pos = Alpha_pos; Delta_pos = Alpha_pos;
for it = 1:MaxIt
    for i = 1:N
        X(i,:) = bound(X(i,:),lb,ub); fit = fobj(X(i,:));
        if fit < Alpha_score
            Delta_score = Beta_score; Delta_pos = Beta_pos;
            Beta_score = Alpha_score; Beta_pos = Alpha_pos;
            Alpha_score = fit; Alpha_pos = X(i,:);
        elseif fit < Beta_score
            Delta_score = Beta_score; Delta_pos = Beta_pos;
            Beta_score = fit; Beta_pos = X(i,:);
        elseif fit < Delta_score
            Delta_score = fit; Delta_pos = X(i,:);
        end
    end
    a = 2 - 2*it/MaxIt;
    for i = 1:N
        X(i,:) = gwoUpdate(X(i,:),Alpha_pos,Beta_pos,Delta_pos,a);
    end
    curve(it) = Alpha_score;
end
end

function Xnew = gwoUpdate(X,Alpha,Beta,Delta,a)
dim = numel(X);
r1=rand(1,dim); r2=rand(1,dim); A1=2*a*r1-a; C1=2*r2; X1=Alpha-A1.*abs(C1.*Alpha-X);
r1=rand(1,dim); r2=rand(1,dim); A2=2*a*r1-a; C2=2*r2; X2=Beta-A2.*abs(C2.*Beta-X);
r1=rand(1,dim); r2=rand(1,dim); A3=2*a*r1-a; C3=2*r2; X3=Delta-A3.*abs(C3.*Delta-X);
Xnew = (X1+X2+X3)/3;
end

function [bestF,bestX,curve] = algWOA(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); bestF = inf; bestX = zeros(1,dim); curve = zeros(MaxIt,1); b = 1;
for it = 1:MaxIt
    for i = 1:N
        X(i,:) = bound(X(i,:),lb,ub); fit = fobj(X(i,:));
        if fit < bestF, bestF = fit; bestX = X(i,:); end
    end
    a = 2 - 2*it/MaxIt;
    for i = 1:N
        A = 2*a*rand-a; C = 2*rand; p = rand; l = -1 + 2*rand;
        if p < 0.5
            if abs(A) < 1
                D = abs(C*bestX-X(i,:)); X(i,:) = bestX-A.*D;
            else
                xr = X(randi(N),:); D = abs(C*xr-X(i,:)); X(i,:) = xr-A.*D;
            end
        else
            D = abs(bestX-X(i,:)); X(i,:) = D.*exp(b*l).*cos(2*pi*l)+bestX;
        end
    end
    curve(it) = bestF;
end
end

function [bestF,bestX,curve] = algDE(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); F = 0.5; CR = 0.9; fit = zeros(N,1); curve = zeros(MaxIt,1);
for i = 1:N, fit(i) = fobj(X(i,:)); end
[bestF,id] = min(fit); bestX = X(id,:);
for it = 1:MaxIt
    for i = 1:N
        idx = randperm(N,3);
        while any(idx==i), idx = randperm(N,3); end
        V = X(idx(1),:) + F*(X(idx(2),:)-X(idx(3),:));
        V = bound(V,lb,ub); U = X(i,:); jrand = randi(dim);
        for j = 1:dim
            if rand < CR || j == jrand, U(j) = V(j); end
        end
        fu = fobj(U);
        if fu < fit(i)
            X(i,:) = U; fit(i) = fu;
            if fu < bestF, bestF = fu; bestX = U; end
        end
    end
    curve(it) = bestF;
end
end

function [bestF,bestX,curve] = algSSA(N,MaxIt,lb,ub,dim,fobj)
X = initRand(N,dim,lb,ub); fit = zeros(N,1); curve = zeros(MaxIt,1); ST = 0.8; PD = 0.2; SD = 0.1;
for i = 1:N, fit(i) = fobj(X(i,:)); end
[bestF,id] = min(fit); bestX = X(id,:); [worstF,wid] = max(fit); worstX = X(wid,:);
for it = 1:MaxIt
    [fit,idx] = sort(fit); X = X(idx,:); bestX = X(1,:); bestF = fit(1); worstX = X(end,:); worstF = fit(end);
    nP = round(PD*N); nS = round(SD*N);
    for i = 1:nP
        if rand < ST, X(i,:) = X(i,:).*exp(-i/(rand*MaxIt+eps)); else, X(i,:) = X(i,:) + randn(1,dim); end
    end
    for i = nP+1:N
        if i > N/2, X(i,:) = randn(1,dim).*exp((worstX-X(i,:))/(i^2)); else, X(i,:) = bestX + abs(X(i,:)-bestX).*randn(1,dim); end
    end
    danger = randperm(N,nS);
    for k = danger
        if fit(k) > bestF
            X(k,:) = bestX + randn(1,dim).*abs(X(k,:)-bestX);
        else
            X(k,:) = X(k,:) + (2*rand(1,dim)-1).*abs(X(k,:)-worstX)/(fit(k)-worstF+eps);
        end
    end
    for i = 1:N
        X(i,:) = bound(X(i,:),lb,ub); fit(i) = fobj(X(i,:));
    end
    [bf,id] = min(fit);
    if bf < bestF, bestF = bf; bestX = X(id,:); end
    curve(it) = bestF;
end
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
                if fobj(bound(Y,lb,ub)) < fobj(X(i,:))
                    X(i,:) = Y;
                elseif fobj(bound(Z,lb,ub)) < fobj(X(i,:))
                    X(i,:) = Z;
                end
            end
        end
    end
    curve(it) = bestF;
end
end
