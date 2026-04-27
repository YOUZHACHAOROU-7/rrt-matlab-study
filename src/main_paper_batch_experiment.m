%% main_paper_batch_experiment.m
% Batch experiment for paper tables and figures.
% Run from the repository root or from the src folder.
% Outputs are saved into ./results_batch

clear; clc; close all; rng(2026);

outDir = fullfile(pwd,'results_batch');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% 1. Joint-space waypoints: q0 -> q1 -> q2 -> q3, unit: rad
qWay = deg2rad([
     0  -45   35    0   25    0;
    25  -25   45   20   10   15;
    45  -10   20   35  -15   25;
    70   20  -10   45  -35   40]);

qLim = deg2rad(repmat([-170 170],6,1));
vMax = deg2rad([90 90 90 120 120 180]);
aMax = deg2rad([180 180 180 240 240 360]);
dh = [0 pi/2 0.1625 0; -0.425 0 0 0; -0.3922 0 0 0; 0 pi/2 0.1333 0; 0 -pi/2 0.0997 0; 0 0 0.0996 0];

%% 2. Optimizer settings
lb = [0.25 0.25 0.25];
ub = [4.50 4.50 4.50];
opt.nPop = 40;
opt.maxIter = 150;
opt.wMax = 0.90;
opt.wMin = 0.35;
opt.c1Max = 2.5;
opt.c1Min = 0.5;
opt.c2Max = 2.5;
opt.c2Min = 0.5;
opt.penalty = 1e4;
opt.levyProb = 0.30;
opt.levyBeta = 1.5;
obj = @(x) fitnessTime(x,qWay,qLim,vMax,aMax,opt.penalty);

nRun = 30;
algNames = {'Fixed-353';'PSO';'Tent-Levy-PSO'};
fixedT = [2.5 2.5 2.5];
fixedF = obj(fixedT);

bestX = zeros(nRun,3,2);       % PSO, IPSO
bestF = zeros(nRun,2);
bestTime = zeros(nRun,2);
finalCurves = zeros(opt.maxIter,2,nRun);
hitIter = zeros(nRun,2);

fprintf('Batch experiment starts: %d repeated runs.\n', nRun);
for r = 1:nRun
    rng(202600+r);
    [x1,f1,c1] = basicPSO(obj,lb,ub,opt);
    [x2,f2,c2] = improvedPSO(obj,lb,ub,opt);
    bestX(r,:,1)=x1; bestX(r,:,2)=x2;
    bestF(r,1)=f1; bestF(r,2)=f2;
    bestTime(r,1)=sum(x1); bestTime(r,2)=sum(x2);
    finalCurves(:,1,r)=c1; finalCurves(:,2,r)=c2;
    hitIter(r,1)=firstHit(c1,f1*1.001);
    hitIter(r,2)=firstHit(c2,f2*1.001);
    fprintf('Run %02d/%02d | PSO %.4f | TLPSO %.4f\n', r,nRun,bestTime(r,1),bestTime(r,2));
end

%% 3. Summary tables
FixedRow = table({'Fixed-353'},fixedT(1),fixedT(2),fixedT(3),sum(fixedT),fixedF,NaN,NaN,NaN, ...
    'VariableNames',{'Algorithm','Best_t1','Best_t2','Best_t3','BestTime','BestFitness','MeanTime','StdTime','MeanHitIter'});

[~,psoBestIdx] = min(bestF(:,1));
[~,ipsoBestIdx] = min(bestF(:,2));
PsoRow = makeRow('PSO',bestX(psoBestIdx,:,1),bestF(psoBestIdx,1),bestTime(:,1),hitIter(:,1));
IpsoRow = makeRow('Tent-Levy-PSO',bestX(ipsoBestIdx,:,2),bestF(ipsoBestIdx,2),bestTime(:,2),hitIter(:,2));
Summary = [FixedRow; PsoRow; IpsoRow];

disp('========== Paper Summary Table ==========');
disp(Summary);
writetable(Summary,fullfile(outDir,'paper_summary_table.csv'));

RunID = (1:nRun)';
Detail = table(RunID,bestTime(:,1),bestF(:,1),hitIter(:,1),bestTime(:,2),bestF(:,2),hitIter(:,2), ...
    'VariableNames',{'RunID','PSO_Time','PSO_Fitness','PSO_HitIter','TLPSO_Time','TLPSO_Fitness','TLPSO_HitIter'});
writetable(Detail,fullfile(outDir,'repeated_runs_detail.csv'));

improveVsFixed = 100*(sum(fixedT)-Summary.BestTime)./sum(fixedT);
ImproveTable = table(Summary.Algorithm,Summary.BestTime,improveVsFixed, ...
    'VariableNames',{'Algorithm','BestTime','ReductionComparedWithFixedPercent'});
writetable(ImproveTable,fullfile(outDir,'improvement_vs_fixed.csv'));

disp('========== Time Reduction Compared With Fixed 3-5-3 ==========' );
disp(ImproveTable);

%% 4. Figures
meanPSO = mean(finalCurves(:,1,:),3);
meanIPSO = mean(finalCurves(:,2,:),3);
fig1=figure('Color','w','Name','Average convergence');
plot(meanPSO,'--','LineWidth',1.6); hold on;
plot(meanIPSO,'-','LineWidth',2.0); grid on;
xlabel('Iteration'); ylabel('Mean best fitness');
legend('PSO','Tent-Levy-PSO','Location','northeast');
title('Average convergence curves over repeated runs');
saveas(fig1,fullfile(outDir,'fig_average_convergence.png'));

% Use the best IPSO trajectory for kinematic curves.
bestTraj = buildTraj353(qWay,bestX(ipsoBestIdx,:,2),260);
fig2=figure('Color','w','Name','Joint displacement');
plot(bestTraj.time,rad2deg(bestTraj.q),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint angle / deg'); title('Joint displacement curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig2,fullfile(outDir,'fig_joint_displacement.png'));

fig3=figure('Color','w','Name','Joint velocity');
plot(bestTraj.time,rad2deg(bestTraj.qd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint velocity / deg/s'); title('Joint velocity curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig3,fullfile(outDir,'fig_joint_velocity.png'));

fig4=figure('Color','w','Name','Joint acceleration');
plot(bestTraj.time,rad2deg(bestTraj.qdd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint acceleration / deg/s^2'); title('Joint acceleration curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig4,fullfile(outDir,'fig_joint_acceleration.png'));

P=zeros(size(bestTraj.q,1),3);
for i=1:size(bestTraj.q,1)
    T=fkineDH(dh,bestTraj.q(i,:));
    P(i,:)=T(1:3,4)';
end
fig5=figure('Color','w','Name','End-effector trajectory');
plot3(P(:,1),P(:,2),P(:,3),'LineWidth',2); grid on; axis equal;
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m'); title('End-effector trajectory');
saveas(fig5,fullfile(outDir,'fig_end_effector_trajectory.png'));

fprintf('\nAll tables and figures have been saved to: %s\n', outDir);

%% ================= local functions =================
function row = makeRow(name,x,f,timeVec,hitVec)
    [bestTime,~] = min(timeVec);
    row = table({name},x(1),x(2),x(3),bestTime,f,mean(timeVec),std(timeVec),mean(hitVec), ...
        'VariableNames',{'Algorithm','Best_t1','Best_t2','Best_t3','BestTime','BestFitness','MeanTime','StdTime','MeanHitIter'});
end

function idx = firstHit(curve,threshold)
    idx = find(curve <= threshold,1,'first');
    if isempty(idx), idx = numel(curve); end
end

function f = fitnessTime(t,qWay,qLim,vMax,aMax,penalty)
    if any(t<=0) || any(~isfinite(t)), f=1e12; return; end
    tr=buildTraj353(qWay,t,150);
    vViol=max(0,abs(tr.qd)-vMax);
    aViol=max(0,abs(tr.qdd)-aMax);
    qViol=max(0,qLim(:,1)'-tr.q)+max(0,tr.q-qLim(:,2)');
    smooth=sum(abs(diff(tr.qdd,1,1)),'all')*1e-4;
    f=sum(t)+penalty*(sum(vViol.^2,'all')+sum(aViol.^2,'all')+sum(qViol.^2,'all'))+smooth;
end

function tr=buildTraj353(qWay,tSeg,nEach)
    nJ=size(qWay,2); q=[]; qd=[]; qdd=[]; time=[]; offset=0;
    for s=1:3
        tt=linspace(0,tSeg(s),nEach)';
        if s>1, tt=tt(2:end); end
        qs=zeros(numel(tt),nJ); qds=qs; qdds=qs;
        for j=1:nJ
            coef=coef353(qWay(:,j),tSeg);
            if s==1, c=coef(1:4); elseif s==2, c=coef(5:10); else, c=coef(11:14); end
            dc=polyder(c); ddc=polyder(dc);
            qs(:,j)=polyval(c,tt); qds(:,j)=polyval(dc,tt); qdds(:,j)=polyval(ddc,tt);
        end
        q=[q;qs]; qd=[qd;qds]; qdd=[qdd;qdds]; time=[time;offset+tt]; offset=offset+tSeg(s);
    end
    tr.q=q; tr.qd=qd; tr.qdd=qdd; tr.time=time;
end

function c=coef353(theta,t)
    t1=t(1); t2=t(2); t3=t(3); th0=theta(1); th1=theta(2); th2=theta(3); th3=theta(4);
    A=zeros(14,14); b=zeros(14,1); r=0;
    r=r+1; A(r,1:4)=basis(3,0,0); b(r)=th0;
    r=r+1; A(r,1:4)=basis(3,0,1);
    r=r+1; A(r,1:4)=basis(3,0,2);
    r=r+1; A(r,1:4)=basis(3,t1,0); b(r)=th1;
    r=r+1; A(r,5:10)=basis(5,0,0); b(r)=th1;
    r=r+1; A(r,1:4)=basis(3,t1,1); A(r,5:10)=-basis(5,0,1);
    r=r+1; A(r,1:4)=basis(3,t1,2); A(r,5:10)=-basis(5,0,2);
    r=r+1; A(r,5:10)=basis(5,t2,0); b(r)=th2;
    r=r+1; A(r,11:14)=basis(3,0,0); b(r)=th2;
    r=r+1; A(r,5:10)=basis(5,t2,1); A(r,11:14)=-basis(3,0,1);
    r=r+1; A(r,5:10)=basis(5,t2,2); A(r,11:14)=-basis(3,0,2);
    r=r+1; A(r,11:14)=basis(3,t3,0); b(r)=th3;
    r=r+1; A(r,11:14)=basis(3,t3,1);
    r=r+1; A(r,11:14)=basis(3,t3,2);
    c=(A\b)';
end

function v=basis(order,t,der)
    p=order:-1:0; v=zeros(1,order+1);
    for k=1:numel(p)
        if p(k)>=der
            m=1; for d=0:der-1, m=m*(p(k)-d); end
            v(k)=m*t^(p(k)-der);
        end
    end
end

function [gb,gf,curve]=basicPSO(obj,lb,ub,opt)
    dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; V=zeros(n,dim); P=X;
    Pf=zeros(n,1); for i=1:n, Pf(i)=obj(X(i,:)); end
    [gf,id]=min(Pf); gb=P(id,:); curve=zeros(opt.maxIter,1);
    for it=1:opt.maxIter
        w=opt.wMax-(opt.wMax-opt.wMin)*it/opt.maxIter;
        for i=1:n
            V(i,:)=w*V(i,:)+2*rand(1,dim).*(P(i,:)-X(i,:))+2*rand(1,dim).*(gb-X(i,:));
            X(i,:)=min(max(X(i,:)+V(i,:),lb),ub); f=obj(X(i,:));
            if f<Pf(i), Pf(i)=f; P(i,:)=X(i,:); end
            if f<gf, gf=f; gb=X(i,:); end
        end
        curve(it)=gf;
    end
end

function [gb,gf,curve]=improvedPSO(obj,lb,ub,opt)
    dim=numel(lb); n=opt.nPop; X=tentInit(n,dim,lb,ub); V=zeros(n,dim); P=X;
    Pf=zeros(n,1); for i=1:n, Pf(i)=obj(X(i,:)); end
    [gf,id]=min(Pf); gb=P(id,:); curve=zeros(opt.maxIter,1); stall=0; last=gf;
    for it=1:opt.maxIter
        w=opt.wMin+(opt.wMax-opt.wMin)*cos(pi*it/(2*opt.maxIter))^2;
        c1=opt.c1Max-(opt.c1Max-opt.c1Min)*it/opt.maxIter;
        c2=opt.c2Min+(opt.c2Max-opt.c2Min)*it/opt.maxIter;
        for i=1:n
            V(i,:)=w*V(i,:)+c1*rand(1,dim).*(P(i,:)-X(i,:))+c2*rand(1,dim).*(gb-X(i,:));
            Xn=X(i,:)+V(i,:);
            if rand<opt.levyProb, Xn=Xn+0.08*levy(dim,opt.levyBeta).*(X(i,:)-gb); end
            Xn=min(max(Xn,lb),ub); f=obj(Xn); X(i,:)=Xn;
            if f<Pf(i), Pf(i)=f; P(i,:)=Xn; end
            if f<gf, gf=f; gb=Xn; end
        end
        if abs(last-gf)<1e-9, stall=stall+1; else, stall=0; last=gf; end
        if stall>20
            [~,ord]=sort(Pf,'descend'); k=max(2,round(0.15*n)); idx=ord(1:k);
            X(idx,:)=tentInit(k,dim,lb,ub);
            for m=1:k
                ii=idx(m); P(ii,:)=X(ii,:); Pf(ii)=obj(X(ii,:));
            end
            [tmp,id]=min(Pf); if tmp<gf, gf=tmp; gb=P(id,:); end
            stall=0;
        end
        curve(it)=gf;
    end
end

function X=tentInit(n,dim,lb,ub)
    z=rand(n,dim); mu=0.499;
    for k=1:8
        idx=z<mu; z(idx)=z(idx)/mu; z(~idx)=(1-z(~idx))/(1-mu);
    end
    X=z.*(ub-lb)+lb;
end

function s=levy(dim,beta)
    sig=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    s=randn(1,dim)*sig./(abs(randn(1,dim)).^(1/beta)+eps);
end

function T=fkineDH(dh,q)
    T=eye(4);
    for i=1:size(dh,1)
        a=dh(i,1); al=dh(i,2); d=dh(i,3); th=q(i)+dh(i,4);
        A=[cos(th) -sin(th)*cos(al) sin(th)*sin(al) a*cos(th); sin(th) cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); 0 sin(al) cos(al) d; 0 0 0 1];
        T=T*A;
    end
end
