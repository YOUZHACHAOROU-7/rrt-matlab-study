%% main_paper_igwo_er7_experiment.m
% ER7-900 paper-level experiment: 3-5-3 interpolation + IGWO time optimization.
% Algorithms: Fixed-353, PSO, GWO, WOA, IGWO.
% Run from repository root: RUN_IGWO_ER7_EXPERIMENT

clear; clc; close all; rng(20260428);

rootDir = pwd;
if endsWith(rootDir,[filesep 'src']), rootDir = fileparts(rootDir); end
outDir = fullfile(rootDir,'results_igwo_er7');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% 1. ER7-900 parameters from manual/product sheet
% Joint ranges, deg: J1 ±170, J2 +100~-135, J3 +200~-75, J4 ±190, J5 ±120, J6 ±360
qLimDeg = [-170 170; -135 100; -75 200; -190 190; -120 120; -360 360];
% Max single-axis velocities, deg/s: [300,255,320,450,450,720]
vMaxDeg = [300 255 320 450 450 720];
% Conservative acceleration limits for simulation, deg/s^2. Can be adjusted if official acceleration is available.
aMaxDeg = 2.2 * vMaxDeg;
qLim = deg2rad(qLimDeg); vMax = deg2rad(vMaxDeg); aMax = deg2rad(aMaxDeg);

% Approximate DH model reconstructed from the ER7-900 workspace/main dimensions drawing.
% This DH model is used for end-effector trajectory visualization, while time optimization is performed in joint space.
dh = [0       pi/2   0.3377  0;
      0.4275  0      0       0;
      0.3760  0      0       0;
      0       pi/2   0.1336  0;
      0      -pi/2   0.0890  0;
      0       0      0.1164  0];

%% 2. Four typical joint-space waypoints for a sorting/handling task, deg
% The values are inside the ER7-900 joint ranges.
qWayDeg = [
      0   -55    45     0    35     0;
     35   -35    70    45    20    60;
     75   -10    35    90   -30   120;
    115    25     5   130   -70   180];
qWay = deg2rad(qWayDeg);

%% 3. Optimization settings
lb = [0.25 0.25 0.25];
ub = [4.50 4.50 4.50];
fixedT = [2.5 2.5 2.5];

opt.nPop = 45; opt.maxIter = 180;
opt.wMax = 0.90; opt.wMin = 0.35;
opt.penalty = 1e4;
opt.levyBeta = 1.5;
opt.cauchyScale = 0.06;
opt.restartRatio = 0.15;

obj = @(t) fitnessTime3(t,qWay,qLim,vMax,aMax,opt.penalty);
fixedF = obj(fixedT);

algNames = {'PSO','GWO','WOA','IGWO'};
nAlg = numel(algNames); nRun = 30;
bestX = zeros(nRun,3,nAlg); bestF = zeros(nRun,nAlg); bestTime = zeros(nRun,nAlg);
curves = zeros(opt.maxIter,nAlg,nRun); hitIter = zeros(nRun,nAlg);

fprintf('ER7-900 IGWO experiment starts: %d repeated runs.\n',nRun);
for r = 1:nRun
    rng(8200+r);
    [x1,f1,c1] = basicPSO(obj,lb,ub,opt);
    [x2,f2,c2] = gwoBasic(obj,lb,ub,opt);
    [x3,f3,c3] = woaBasic(obj,lb,ub,opt);
    [x4,f4,c4] = igwoPaper(obj,lb,ub,opt);
    Xs = [x1; x2; x3; x4]; Fs = [f1 f2 f3 f4]; Cs = {c1,c2,c3,c4};
    for a = 1:nAlg
        bestX(r,:,a) = Xs(a,:); bestF(r,a) = Fs(a); bestTime(r,a) = sum(Xs(a,:)); curves(:,a,r) = Cs{a};
        hitIter(r,a) = firstHit(Cs{a}, Fs(a)*1.001);
    end
    fprintf('Run %02d/%02d | PSO %.4f | GWO %.4f | WOA %.4f | IGWO %.4f\n', ...
        r,nRun,bestTime(r,1),bestTime(r,2),bestTime(r,3),bestTime(r,4));
end

%% 4. Summary tables
FixedRow = table({'Fixed-353'},fixedT(1),fixedT(2),fixedT(3),sum(fixedT),fixedF,NaN,NaN,NaN, ...
    'VariableNames',{'Algorithm','Best_t1','Best_t2','Best_t3','BestTime','BestFitness','MeanTime','StdTime','MeanHitIter'});
Summary = FixedRow;
for a = 1:nAlg
    [bf,idx] = min(bestF(:,a)); bx = bestX(idx,:,a); bt = sum(bx);
    row = table(algNames(a),bx(1),bx(2),bx(3),bt,bf,mean(bestTime(:,a)),std(bestTime(:,a)),mean(hitIter(:,a)), ...
        'VariableNames',{'Algorithm','Best_t1','Best_t2','Best_t3','BestTime','BestFitness','MeanTime','StdTime','MeanHitIter'});
    Summary = [Summary; row]; %#ok<AGROW>
end

disp('========== ER7-900 IGWO Paper Summary Table =========='); disp(Summary);
writetable(Summary,fullfile(outDir,'paper_summary_igwo_er7.csv'));

ReductionVsFixed = 100*(sum(fixedT)-Summary.BestTime)/sum(fixedT);
ReductionTable = table(Summary.Algorithm,Summary.BestTime,ReductionVsFixed, ...
    'VariableNames',{'Algorithm','BestTime','ReductionComparedWithFixedPercent'});
disp('========== Reduction Compared With Fixed 3-5-3 =========='); disp(ReductionTable);
writetable(ReductionTable,fullfile(outDir,'reduction_vs_fixed_igwo_er7.csv'));

RunID = (1:nRun)';
Detail = table(RunID,bestTime(:,1),bestF(:,1),bestTime(:,2),bestF(:,2),bestTime(:,3),bestF(:,3),bestTime(:,4),bestF(:,4), ...
    'VariableNames',{'RunID','PSO_Time','PSO_Fit','GWO_Time','GWO_Fit','WOA_Time','WOA_Fit','IGWO_Time','IGWO_Fit'});
writetable(Detail,fullfile(outDir,'repeated_runs_igwo_er7.csv'));

%% 5. Ablation experiment for IGWO strategies
fprintf('\nAblation experiment starts.\n');
abNames = {'GWO','Cos-GWO','Cos-Mem-GWO','IGWO'}; nAb = numel(abNames); abRun = 20;
abTime = zeros(abRun,nAb); abFit = zeros(abRun,nAb); abCurves = zeros(opt.maxIter,nAb,abRun);
for r = 1:abRun
    rng(9200+r);
    [x0,f0,c0] = gwoBasic(obj,lb,ub,opt);
    [x1,f1,c1] = igwoVariant(obj,lb,ub,opt,1);
    [x2,f2,c2] = igwoVariant(obj,lb,ub,opt,2);
    [x3,f3,c3] = igwoVariant(obj,lb,ub,opt,3);
    Xs=[x0;x1;x2;x3]; Fs=[f0 f1 f2 f3]; Cs={c0,c1,c2,c3};
    for k=1:nAb
        abTime(r,k)=sum(Xs(k,:)); abFit(r,k)=Fs(k); abCurves(:,k,r)=Cs{k};
    end
    fprintf('Ablation %02d/%02d | GWO %.4f | Cos %.4f | CosMem %.4f | IGWO %.4f\n',r,abRun,abTime(r,1),abTime(r,2),abTime(r,3),abTime(r,4));
end
Ablation = table(abNames',min(abTime)',mean(abTime)',std(abTime)',mean(abFit)', ...
    'VariableNames',{'Algorithm','BestTime','MeanTime','StdTime','MeanFitness'});
disp('========== IGWO Ablation Table =========='); disp(Ablation);
writetable(Ablation,fullfile(outDir,'ablation_igwo_er7.csv'));

%% 6. Figures
fig1=figure('Color','w','Name','ER7 average convergence'); hold on;
styles={'--',':','-.','-'};
for a=1:nAlg, plot(mean(curves(:,a,:),3),styles{a},'LineWidth',1.8); end
grid on; xlabel('Iteration'); ylabel('Mean best fitness'); legend(algNames,'Location','northeast'); title('Average convergence curves');
saveas(fig1,fullfile(outDir,'fig_er7_average_convergence.png'));

figAb=figure('Color','w','Name','IGWO ablation convergence'); hold on;
for a=1:nAb, plot(mean(abCurves(:,a,:),3),'LineWidth',1.6); end
grid on; xlabel('Iteration'); ylabel('Mean best fitness'); legend(abNames,'Location','northeast'); title('Ablation convergence curves');
saveas(figAb,fullfile(outDir,'fig_er7_ablation_convergence.png'));

[~,idxIGWO] = min(bestF(:,4)); bestT = bestX(idxIGWO,:,4); traj = buildTraj353(qWay,bestT,260);

fig2=figure('Color','w','Name','ER7 joint displacement'); plot(traj.time,rad2deg(traj.q),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint angle / deg'); title('Joint displacement curves'); legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig2,fullfile(outDir,'fig_er7_joint_displacement.png'));

fig3=figure('Color','w','Name','ER7 joint velocity'); plot(traj.time,rad2deg(traj.qd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint velocity / deg/s'); title('Joint velocity curves'); legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig3,fullfile(outDir,'fig_er7_joint_velocity.png'));

fig4=figure('Color','w','Name','ER7 joint acceleration'); plot(traj.time,rad2deg(traj.qdd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint acceleration / deg/s^2'); title('Joint acceleration curves'); legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig4,fullfile(outDir,'fig_er7_joint_acceleration.png'));

P=zeros(size(traj.q,1),3); for i=1:size(traj.q,1), T=fkineDH(dh,traj.q(i,:)); P(i,:)=T(1:3,4)'; end
fig5=figure('Color','w','Name','ER7 end-effector trajectory'); plot3(P(:,1),P(:,2),P(:,3),'LineWidth',2); grid on; axis equal;
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m'); title('End-effector trajectory');
saveas(fig5,fullfile(outDir,'fig_er7_end_effector_trajectory.png'));

save(fullfile(outDir,'workspace_igwo_er7.mat'));
fprintf('\nAll ER7-900 IGWO tables and figures have been saved to: %s\n',outDir);

%% ================= local functions =================
function idx = firstHit(curve,threshold)
idx=find(curve<=threshold,1,'first'); if isempty(idx), idx=numel(curve); end
end
function f = fitnessTime3(t,qWay,qLim,vMax,aMax,penalty)
if any(t<=0)||any(~isfinite(t)), f=1e12; return; end
tr=buildTraj353(qWay,t,160);
vViol=max(0,abs(tr.qd)-vMax); aViol=max(0,abs(tr.qdd)-aMax);
qViol=max(0,qLim(:,1)'-tr.q)+max(0,tr.q-qLim(:,2)');
smooth=sum(abs(diff(tr.qdd,1,1)),'all')*1e-4;
f=sum(t)+penalty*(sum(vViol.^2,'all')+sum(aViol.^2,'all')+sum(qViol.^2,'all'))+smooth;
end
function tr=buildTraj353(qWay,tSeg,nEach)
nJ=size(qWay,2); q=[]; qd=[]; qdd=[]; time=[]; offset=0;
for s=1:3
    tt=linspace(0,tSeg(s),nEach)'; if s>1, tt=tt(2:end); end
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
r=r+1; A(r,1:4)=basis(3,0,0); b(r)=th0; r=r+1; A(r,1:4)=basis(3,0,1); r=r+1; A(r,1:4)=basis(3,0,2);
r=r+1; A(r,1:4)=basis(3,t1,0); b(r)=th1; r=r+1; A(r,5:10)=basis(5,0,0); b(r)=th1;
r=r+1; A(r,1:4)=basis(3,t1,1); A(r,5:10)=-basis(5,0,1); r=r+1; A(r,1:4)=basis(3,t1,2); A(r,5:10)=-basis(5,0,2);
r=r+1; A(r,5:10)=basis(5,t2,0); b(r)=th2; r=r+1; A(r,11:14)=basis(3,0,0); b(r)=th2;
r=r+1; A(r,5:10)=basis(5,t2,1); A(r,11:14)=-basis(3,0,1); r=r+1; A(r,5:10)=basis(5,t2,2); A(r,11:14)=-basis(3,0,2);
r=r+1; A(r,11:14)=basis(3,t3,0); b(r)=th3; r=r+1; A(r,11:14)=basis(3,t3,1); r=r+1; A(r,11:14)=basis(3,t3,2);
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
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; V=zeros(n,dim); P=X; Pf=zeros(n,1);
for i=1:n, Pf(i)=obj(X(i,:)); end
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
function [Alpha_pos,Alpha_score,curve]=gwoBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; [Alpha_pos,Alpha_score,curve]=gwoLoop(obj,lb,ub,opt,X,0);
end
function [bestX,bestF,curve]=woaBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; bestX=zeros(1,dim); bestF=inf; curve=zeros(opt.maxIter,1); b=1;
for it=1:opt.maxIter
    for i=1:n, X(i,:)=min(max(X(i,:),lb),ub); fit=obj(X(i,:)); if fit<bestF, bestF=fit; bestX=X(i,:); end, end
    a=2-2*it/opt.maxIter;
    for i=1:n
        A=2*a*rand-a; C=2*rand; p=rand; l=-1+2*rand;
        if p<0.5
            if abs(A)<1, D=abs(C*bestX-X(i,:)); X(i,:)=bestX-A.*D; else, xr=X(randi(n),:); D=abs(C*xr-X(i,:)); X(i,:)=xr-A.*D; end
        else
            D=abs(bestX-X(i,:)); X(i,:)=D.*exp(b*l).*cos(2*pi*l)+bestX;
        end
    end
    curve(it)=bestF;
end
end
function [best,bf,curve]=igwoPaper(obj,lb,ub,opt)
[best,bf,curve]=igwoVariant(obj,lb,ub,opt,3);
end
function [Alpha_pos,Alpha_score,curve]=igwoVariant(obj,lb,ub,opt,mode)
dim=numel(lb); n=opt.nPop; X=tentInit(n,dim,lb,ub); P=X; Pf=inf(n,1); V=zeros(n,dim);
for i=1:n, Pf(i)=obj(X(i,:)); end
[Alpha_score,id]=min(Pf); Alpha_pos=P(id,:); curve=zeros(opt.maxIter,1);
Beta_pos=Alpha_pos; Delta_pos=Alpha_pos; Beta_score=inf; Delta_score=inf;
for it=1:opt.maxIter
    for i=1:n
        X(i,:)=min(max(X(i,:),lb),ub); fit=obj(X(i,:));
        if fit<Pf(i), Pf(i)=fit; P(i,:)=X(i,:); end
        if fit<Alpha_score
            Delta_score=Beta_score; Delta_pos=Beta_pos; Beta_score=Alpha_score; Beta_pos=Alpha_pos; Alpha_score=fit; Alpha_pos=X(i,:);
        elseif fit<Beta_score
            Delta_score=Beta_score; Delta_pos=Beta_pos; Beta_score=fit; Beta_pos=X(i,:);
        elseif fit<Delta_score
            Delta_score=fit; Delta_pos=X(i,:);
        end
    end
    if mode>=1, a=2*cos(pi*it/(2*opt.maxIter))^2; else, a=2-2*it/opt.maxIter; end
    for i=1:n
        Xg = gwoPosition(X(i,:),Alpha_pos,Beta_pos,Delta_pos,a);
        if mode>=2
            w=opt.wMin+(opt.wMax-opt.wMin)*cos(pi*it/(2*opt.maxIter))^2;
            V(i,:)=w*V(i,:)+1.2*rand(1,dim).*(P(i,:)-X(i,:))+1.6*rand(1,dim).*(Alpha_pos-X(i,:));
            Xg=0.70*Xg+0.30*(X(i,:)+V(i,:));
        end
        Xg=min(max(Xg,lb),ub);
        if mode>=3 && rand<0.35
            Xu=Alpha_pos+opt.cauchyScale*cauchyRand(1,dim).*(ub-lb).*(1-it/opt.maxIter);
            Xu=min(max(Xu,lb),ub); fu=obj(Xu);
            if fu<Alpha_score, Alpha_score=fu; Alpha_pos=Xu; end
        end
        X(i,:)=Xg;
    end
    curve(it)=Alpha_score;
end
end
function Xnew=gwoPosition(X,Alpha,Beta,Delta,a)
dim=numel(X); r1=rand(1,dim); r2=rand(1,dim); A1=2*a*r1-a; C1=2*r2; X1=Alpha-A1.*abs(C1.*Alpha-X);
r1=rand(1,dim); r2=rand(1,dim); A2=2*a*r1-a; C2=2*r2; X2=Beta-A2.*abs(C2.*Beta-X);
r1=rand(1,dim); r2=rand(1,dim); A3=2*a*r1-a; C3=2*r2; X3=Delta-A3.*abs(C3.*Delta-X);
Xnew=(X1+X2+X3)/3;
end
function X=tentInit(n,dim,lb,ub)
z=rand(n,dim); mu=0.499;
for k=1:8, idx=z<mu; z(idx)=z(idx)/mu; z(~idx)=(1-z(~idx))/(1-mu); end
X=z.*(ub-lb)+lb;
end
function r=cauchyRand(m,n), r=tan(pi*(rand(m,n)-0.5)); r=max(min(r,10),-10); end
function T=fkineDH(dh,q)
T=eye(4);
for i=1:size(dh,1)
    a=dh(i,1); al=dh(i,2); d=dh(i,3); th=q(i)+dh(i,4);
    A=[cos(th) -sin(th)*cos(al) sin(th)*sin(al) a*cos(th); sin(th) cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); 0 sin(al) cos(al) d; 0 0 0 1];
    T=T*A;
end
end
function [Alpha_pos,Alpha_score,curve]=gwoLoop(obj,lb,ub,opt,X,unused)
dim=numel(lb); n=size(X,1); Alpha_score=inf; Beta_score=inf; Delta_score=inf; Alpha_pos=zeros(1,dim); Beta_pos=Alpha_pos; Delta_pos=Alpha_pos; curve=zeros(opt.maxIter,1);
for it=1:opt.maxIter
    for i=1:n
        X(i,:)=min(max(X(i,:),lb),ub); fit=obj(X(i,:));
        if fit<Alpha_score
            Delta_score=Beta_score; Delta_pos=Beta_pos; Beta_score=Alpha_score; Beta_pos=Alpha_pos; Alpha_score=fit; Alpha_pos=X(i,:);
        elseif fit<Beta_score
            Delta_score=Beta_score; Delta_pos=Beta_pos; Beta_score=fit; Beta_pos=X(i,:);
        elseif fit<Delta_score
            Delta_score=fit; Delta_pos=X(i,:);
        end
    end
    a=2-2*it/opt.maxIter;
    for i=1:n, X(i,:)=gwoPosition(X(i,:),Alpha_pos,Beta_pos,Delta_pos,a); end
    curve(it)=Alpha_score;
end
end
