%% main_er7_mo_trajectory_experiment.m
% ER7-900 multi-objective trajectory optimization
% Objective = normalized time + normalized jerk + normalized energy proxy
% Algorithms: Fixed-353, PSO, GWO, WOA, SSA, HHO, DE, IGWO

clear; clc; close all; rng(20260509);
rootDir = pwd;
if endsWith(rootDir,[filesep 'src']), rootDir = fileparts(rootDir); end
outDir = fullfile(rootDir,'results_er7_multiobjective');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ER7-900 constraints
qLimDeg = [-170 170; -135 100; -75 200; -190 190; -120 120; -360 360];
vMaxDeg = [300 255 320 450 450 720];
aMaxDeg = 2.2 * vMaxDeg;
qLim = deg2rad(qLimDeg); vMax = deg2rad(vMaxDeg); aMax = deg2rad(aMaxDeg);

% Approximate DH model from ER7-900 drawing, used for end-effector visualization only.
dh = [0       pi/2   0.3377  0;
      0.4275  0      0       0;
      0.3760  0      0       0;
      0       pi/2   0.1336  0;
      0      -pi/2   0.0890  0;
      0       0      0.1164  0];

%% Typical lightweight refractory part handling waypoints, deg
qWayDeg = [
      0   -55    45     0    35     0;
     35   -35    70    45    20    60;
     75   -10    35    90   -30   120;
    115    25     5   130   -70   180];
qWay = deg2rad(qWayDeg);

%% Optimization setting
lb = [0.25 0.25 0.25]; ub = [4.50 4.50 4.50]; fixedT = [2.5 2.5 2.5];
opt.nPop = 45; opt.maxIter = 180; opt.penalty = 1e5;
opt.wTime = 0.50; opt.wJerk = 0.25; opt.wEnergy = 0.25;
opt.fixedMetric = trajectoryMetrics(fixedT,qWay,qLim,vMax,aMax);
obj = @(t) multiObjFitness(t,qWay,qLim,vMax,aMax,opt);
fixedFit = obj(fixedT); fixedM = trajectoryMetrics(fixedT,qWay,qLim,vMax,aMax);

algNames = {'PSO','GWO','WOA','SSA','HHO','DE','IGWO'}; nAlg=numel(algNames); nRun=30;
bestX=zeros(nRun,3,nAlg); bestF=zeros(nRun,nAlg); curves=zeros(opt.maxIter,nAlg,nRun);
metricsTime=zeros(nRun,nAlg); metricsJerk=zeros(nRun,nAlg); metricsEnergy=zeros(nRun,nAlg);

fprintf('ER7 multi-objective experiment starts: %d repeated runs.\n',nRun);
for r=1:nRun
    rng(12000+r);
    [x1,f1,c1]=basicPSO(obj,lb,ub,opt);
    [x2,f2,c2]=gwoBasic(obj,lb,ub,opt);
    [x3,f3,c3]=woaBasic(obj,lb,ub,opt);
    [x4,f4,c4]=ssaBasic(obj,lb,ub,opt);
    [x5,f5,c5]=hhoBasic(obj,lb,ub,opt);
    [x6,f6,c6]=deBasic(obj,lb,ub,opt);
    [x7,f7,c7]=igwo3(obj,lb,ub,opt);
    Xs=[x1;x2;x3;x4;x5;x6;x7]; Fs=[f1 f2 f3 f4 f5 f6 f7]; Cs={c1,c2,c3,c4,c5,c6,c7};
    for a=1:nAlg
        bestX(r,:,a)=Xs(a,:); bestF(r,a)=Fs(a); curves(:,a,r)=Cs{a};
        m=trajectoryMetrics(Xs(a,:),qWay,qLim,vMax,aMax);
        metricsTime(r,a)=m.time; metricsJerk(r,a)=m.jerk; metricsEnergy(r,a)=m.energy;
    end
    fprintf('Run %02d/%02d | PSO %.4f | GWO %.4f | WOA %.4f | SSA %.4f | HHO %.4f | DE %.4f | IGWO %.4f\n',r,nRun,Fs);
end

%% Summary
Summary = table({'Fixed-353'},fixedT(1),fixedT(2),fixedT(3),fixedM.time,fixedM.jerk,fixedM.energy,fixedFit,NaN,NaN,NaN, ...
    'VariableNames',{'Algorithm','t1','t2','t3','Time','Jerk','EnergyProxy','BestFitness','MeanFitness','StdFitness','MeanTime'});
for a=1:nAlg
    [bf,idx]=min(bestF(:,a)); bx=bestX(idx,:,a); bm=trajectoryMetrics(bx,qWay,qLim,vMax,aMax);
    row=table(algNames(a),bx(1),bx(2),bx(3),bm.time,bm.jerk,bm.energy,bf,mean(bestF(:,a)),std(bestF(:,a)),mean(metricsTime(:,a)), ...
        'VariableNames',Summary.Properties.VariableNames);
    Summary=[Summary;row]; %#ok<AGROW>
end

disp('========== ER7 Multi-objective Summary Table =========='); disp(Summary);
writetable(Summary,fullfile(outDir,'er7_multiobjective_summary.csv'));

Reduction = table(Summary.Algorithm,Summary.Time,100*(fixedM.time-Summary.Time)/fixedM.time,Summary.Jerk,100*(fixedM.jerk-Summary.Jerk)/fixedM.jerk,Summary.EnergyProxy,100*(fixedM.energy-Summary.EnergyProxy)/fixedM.energy, ...
    'VariableNames',{'Algorithm','Time','TimeReductionPercent','Jerk','JerkReductionPercent','EnergyProxy','EnergyReductionPercent'});
disp('========== Reduction Compared With Fixed 3-5-3 =========='); disp(Reduction);
writetable(Reduction,fullfile(outDir,'er7_multiobjective_reduction.csv'));

%% Ablation
abNames={'GWO','Tent-GWO','Tent-Cos-GWO','IGWO'}; abRun=20; nAb=numel(abNames);
abF=zeros(abRun,nAb); abT=zeros(abRun,nAb); abJ=zeros(abRun,nAb); abE=zeros(abRun,nAb); abCurves=zeros(opt.maxIter,nAb,abRun);
fprintf('\nER7 multi-objective ablation starts.\n');
for r=1:abRun
    rng(15000+r);
    [xa,fa,ca]=gwoBasic(obj,lb,ub,opt);
    [xb,fb,cb]=igwoVariant(obj,lb,ub,opt,1);
    [xc,fc,cc]=igwoVariant(obj,lb,ub,opt,2);
    [xd,fd,cd]=igwoVariant(obj,lb,ub,opt,3);
    Xs=[xa;xb;xc;xd]; Fs=[fa fb fc fd]; Cs={ca,cb,cc,cd};
    for k=1:nAb
        mm=trajectoryMetrics(Xs(k,:),qWay,qLim,vMax,aMax);
        abF(r,k)=Fs(k); abT(r,k)=mm.time; abJ(r,k)=mm.jerk; abE(r,k)=mm.energy; abCurves(:,k,r)=Cs{k};
    end
    fprintf('Ablation %02d/%02d | GWO %.4f | Tent %.4f | TentCos %.4f | IGWO %.4f\n',r,abRun,Fs);
end
Ablation=table(abNames',min(abF)',mean(abF)',std(abF)',mean(abT)',mean(abJ)',mean(abE)', ...
    'VariableNames',{'Algorithm','BestFitness','MeanFitness','StdFitness','MeanTime','MeanJerk','MeanEnergyProxy'});
disp('========== ER7 Multi-objective Ablation Table =========='); disp(Ablation);
writetable(Ablation,fullfile(outDir,'er7_multiobjective_ablation.csv'));

%% Figures
fig=figure('Color','w','Name','ER7 multi-objective convergence'); hold on;
for a=1:nAlg, semilogy(mean(curves(:,a,:),3),'LineWidth',1.5); end
grid on; xlabel('Iteration'); ylabel('Mean best fitness'); title('Average convergence curves'); legend(algNames,'Location','northeastoutside');
saveas(fig,fullfile(outDir,'fig_mo_average_convergence.png'));

figA=figure('Color','w','Name','ER7 multi-objective ablation'); hold on;
for a=1:nAb, semilogy(mean(abCurves(:,a,:),3),'LineWidth',1.5); end
grid on; xlabel('Iteration'); ylabel('Mean best fitness'); title('Ablation convergence curves'); legend(abNames,'Location','northeast');
saveas(figA,fullfile(outDir,'fig_mo_ablation_convergence.png'));

idxIGWO=find(strcmp(Summary.Algorithm,'IGWO')); bestT=[Summary.t1(idxIGWO) Summary.t2(idxIGWO) Summary.t3(idxIGWO)];
traj=buildTraj353(qWay,bestT,300);
fig2=figure('Color','w','Name','MO joint displacement'); plot(traj.time,rad2deg(traj.q),'LineWidth',1.3); grid on; xlabel('Time / s'); ylabel('Joint angle / deg'); title('Joint displacement curves'); legend('q1','q2','q3','q4','q5','q6','Location','best'); saveas(fig2,fullfile(outDir,'fig_mo_joint_displacement.png'));
fig3=figure('Color','w','Name','MO joint velocity'); plot(traj.time,rad2deg(traj.qd),'LineWidth',1.3); grid on; xlabel('Time / s'); ylabel('Joint velocity / deg/s'); title('Joint velocity curves'); legend('q1','q2','q3','q4','q5','q6','Location','best'); saveas(fig3,fullfile(outDir,'fig_mo_joint_velocity.png'));
fig4=figure('Color','w','Name','MO joint acceleration'); plot(traj.time,rad2deg(traj.qdd),'LineWidth',1.3); grid on; xlabel('Time / s'); ylabel('Joint acceleration / deg/s^2'); title('Joint acceleration curves'); legend('q1','q2','q3','q4','q5','q6','Location','best'); saveas(fig4,fullfile(outDir,'fig_mo_joint_acceleration.png'));
fig5=figure('Color','w','Name','MO joint jerk'); plot(traj.time,rad2deg(traj.qddd),'LineWidth',1.3); grid on; xlabel('Time / s'); ylabel('Joint jerk / deg/s^3'); title('Joint jerk curves'); legend('q1','q2','q3','q4','q5','q6','Location','best'); saveas(fig5,fullfile(outDir,'fig_mo_joint_jerk.png'));
P=zeros(size(traj.q,1),3); for i=1:size(traj.q,1), T=fkineDH(dh,traj.q(i,:)); P(i,:)=T(1:3,4)'; end
fig6=figure('Color','w','Name','MO end-effector trajectory'); plot3(P(:,1),P(:,2),P(:,3),'LineWidth',2); grid on; axis equal; xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m'); title('End-effector trajectory'); saveas(fig6,fullfile(outDir,'fig_mo_end_effector_trajectory.png'));

fprintf('\nAll ER7 multi-objective tables and figures have been saved to: %s\n',outDir);

%% ================= functions =================
function F=multiObjFitness(t,qWay,qLim,vMax,aMax,opt)
if any(t<=0)||any(~isfinite(t)), F=1e12; return; end
m=trajectoryMetrics(t,qWay,qLim,vMax,aMax); base=opt.fixedMetric;
tr=buildTraj353(qWay,t,160);
vViol=max(0,abs(tr.qd)-vMax); aViol=max(0,abs(tr.qdd)-aMax); qViol=max(0,qLim(:,1)'-tr.q)+max(0,tr.q-qLim(:,2)');
pen=opt.penalty*(sum(vViol.^2,'all')+sum(aViol.^2,'all')+sum(qViol.^2,'all'));
F=opt.wTime*(m.time/base.time)+opt.wJerk*(m.jerk/base.jerk)+opt.wEnergy*(m.energy/base.energy)+pen;
end
function m=trajectoryMetrics(t,qWay,qLim,vMax,aMax)
tr=buildTraj353(qWay,t,220);
m.time=sum(t);
m.jerk=trapz(tr.time,sum(tr.qddd.^2,2));
m.energy=trapz(tr.time,sum(tr.qd.^2,2)+0.05*sum(tr.qdd.^2,2));
end
function tr=buildTraj353(qWay,tSeg,nEach)
nJ=size(qWay,2); q=[]; qd=[]; qdd=[]; qddd=[]; time=[]; offset=0;
for s=1:3
    tt=linspace(0,tSeg(s),nEach)'; if s>1, tt=tt(2:end); end
    qs=zeros(numel(tt),nJ); qds=qs; qdds=qs; qddds=qs;
    for j=1:nJ
        coef=coef353(qWay(:,j),tSeg);
        if s==1, c=coef(1:4); elseif s==2, c=coef(5:10); else, c=coef(11:14); end
        dc=polyder(c); ddc=polyder(dc); dddc=polyder(ddc);
        qs(:,j)=polyval(c,tt); qds(:,j)=polyval(dc,tt); qdds(:,j)=polyval(ddc,tt); qddds(:,j)=polyval(dddc,tt);
    end
    q=[q;qs]; qd=[qd;qds]; qdd=[qdd;qdds]; qddd=[qddd;qddds]; time=[time;offset+tt]; offset=offset+tSeg(s);
end
tr.q=q; tr.qd=qd; tr.qdd=qdd; tr.qddd=qddd; tr.time=time;
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
p=order:-1:0; v=zeros(1,order+1); for k=1:numel(p), if p(k)>=der, m=1; for d=0:der-1, m=m*(p(k)-d); end, v(k)=m*t^(p(k)-der); end, end
end
function [best,bf,curve]=basicPSO(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; V=zeros(n,dim); P=X; Pf=arrayfun(@(i)obj(X(i,:)),1:n)'; [bf,id]=min(Pf); best=P(id,:); curve=zeros(opt.maxIter,1);
for it=1:opt.maxIter, w=0.9-0.5*it/opt.maxIter; for i=1:n, V(i,:)=w*V(i,:)+2*rand(1,dim).*(P(i,:)-X(i,:))+2*rand(1,dim).*(best-X(i,:)); X(i,:)=bound(X(i,:)+V(i,:),lb,ub); f=obj(X(i,:)); if f<Pf(i), Pf(i)=f; P(i,:)=X(i,:); end, if f<bf, bf=f; best=X(i,:); end, end, curve(it)=bf; end
end
function [best,bf,curve]=gwoBasic(obj,lb,ub,opt)
X=rand(opt.nPop,numel(lb)).*(ub-lb)+lb; [best,bf,curve]=gwoCore(obj,lb,ub,opt,X,0);
end
function [best,bf,curve]=igwo3(obj,lb,ub,opt), [best,bf,curve]=igwoVariant(obj,lb,ub,opt,3); end
function [best,bf,curve]=igwoVariant(obj,lb,ub,opt,mode)
dim=numel(lb); n=opt.nPop; if mode>=1, X=tentInit(n,dim,lb,ub); else, X=rand(n,dim).*(ub-lb)+lb; end
[best,bf,curve]=gwoCore(obj,lb,ub,opt,X,mode);
end
function [Alpha,AlphaScore,curve]=gwoCore(obj,lb,ub,opt,X,mode)
[n,dim]=size(X); Alpha=zeros(1,dim); Beta=Alpha; Delta=Alpha; AlphaScore=inf; BetaScore=inf; DeltaScore=inf; curve=zeros(opt.maxIter,1);
for it=1:opt.maxIter
    for i=1:n
        X(i,:)=bound(X(i,:),lb,ub); fit=obj(X(i,:));
        if fit<AlphaScore, DeltaScore=BetaScore; Delta=Beta; BetaScore=AlphaScore; Beta=Alpha; AlphaScore=fit; Alpha=X(i,:);
        elseif fit<BetaScore, DeltaScore=BetaScore; Delta=Beta; BetaScore=fit; Beta=X(i,:);
        elseif fit<DeltaScore, DeltaScore=fit; Delta=X(i,:); end
    end
    if mode>=2, a=2*cos((pi/2)*(it/opt.maxIter)); else, a=2-2*it/opt.maxIter; end
    for i=1:n
        Xn=gwoPosition(X(i,:),Alpha,Beta,Delta,a); Xn=bound(Xn,lb,ub);
        if mode>=3 && rand < 0.35*(1-it/opt.maxIter)+0.05
            LF=levy(dim,1.5); CY=tan(pi*(rand(1,dim)-0.5)); CY=max(min(CY,10),-10);
            Xm=Xn+0.15*(1-it/opt.maxIter)*LF.*(Alpha-Xn)+0.01*(1-it/opt.maxIter)*CY.*(ub-lb); Xm=bound(Xm,lb,ub);
            if obj(Xm)<obj(Xn), Xn=Xm; end
        end
        X(i,:)=Xn;
    end
    curve(it)=AlphaScore;
end
end
function Xn=gwoPosition(X,A,B,D,a)
dim=numel(X); Xs=zeros(3,dim); leaders=[A;B;D]; for k=1:3, r1=rand(1,dim); r2=rand(1,dim); Ak=2*a*r1-a; Ck=2*r2; Xs(k,:)=leaders(k,:)-Ak.*abs(Ck.*leaders(k,:)-X); end, Xn=mean(Xs,1);
end
function [best,bf,curve]=woaBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; best=X(1,:); bf=inf; curve=zeros(opt.maxIter,1); b=1;
for it=1:opt.maxIter, for i=1:n, X(i,:)=bound(X(i,:),lb,ub); f=obj(X(i,:)); if f<bf, bf=f; best=X(i,:); end, end, a=2-2*it/opt.maxIter; for i=1:n, A=2*a*rand-a; C=2*rand; p=rand; l=-1+2*rand; if p<0.5, if abs(A)<1, X(i,:)=best-A.*abs(C*best-X(i,:)); else, xr=X(randi(n),:); X(i,:)=xr-A.*abs(C*xr-X(i,:)); end, else, X(i,:)=abs(best-X(i,:)).*exp(b*l).*cos(2*pi*l)+best; end, end, curve(it)=bf; end
end
function [best,bf,curve]=ssaBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; curve=zeros(opt.maxIter,1); ST=0.8; PD=0.2;
fit=arrayfun(@(i)obj(X(i,:)),1:n)'; [bf,id]=min(fit); best=X(id,:);
for it=1:opt.maxIter, [fit,idx]=sort(fit); X=X(idx,:); for i=1:n, if i<=round(PD*n), if rand<ST, X(i,:)=X(i,:).*exp(-i/(rand*opt.maxIter+eps)); else, X(i,:)=X(i,:)+randn(1,dim); end, else, X(i,:)=X(i,:)+randn(1,dim).*abs(X(i,:)-best); end, X(i,:)=bound(X(i,:),lb,ub); fit(i)=obj(X(i,:)); end, [tmp,id]=min(fit); if tmp<bf, bf=tmp; best=X(id,:); end, curve(it)=bf; end
end
function [best,bf,curve]=hhoBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; bf=inf; best=X(1,:); curve=zeros(opt.maxIter,1);
for it=1:opt.maxIter, for i=1:n, X(i,:)=bound(X(i,:),lb,ub); f=obj(X(i,:)); if f<bf, bf=f; best=X(i,:); end, end, E1=2*(1-it/opt.maxIter); Xm=mean(X,1); for i=1:n, E=E1*(2*rand-1); if abs(E)>=1, Xrand=X(randi(n),:); X(i,:)=Xrand-rand(1,dim).*abs(Xrand-2*rand(1,dim).*X(i,:)); else, X(i,:)=best-E.*abs(best-X(i,:))+0.01*randn(1,dim).*(Xm-X(i,:)); end, end, curve(it)=bf; end
end
function [best,bf,curve]=deBasic(obj,lb,ub,opt)
dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; F=0.5; CR=0.9; fit=arrayfun(@(i)obj(X(i,:)),1:n)'; [bf,id]=min(fit); best=X(id,:); curve=zeros(opt.maxIter,1);
for it=1:opt.maxIter, for i=1:n, ids=randperm(n,3); V=X(ids(1),:)+F*(X(ids(2),:)-X(ids(3),:)); U=X(i,:); jrand=randi(dim); for j=1:dim, if rand<CR||j==jrand, U(j)=V(j); end, end, U=bound(U,lb,ub); fu=obj(U); if fu<fit(i), X(i,:)=U; fit(i)=fu; if fu<bf, bf=fu; best=U; end, end, end, curve(it)=bf; end
end
function X=tentInit(n,dim,lb,ub)
z=rand(n,dim); for k=1:8, idx=z<0.5; z(idx)=z(idx)/0.5; z(~idx)=(1-z(~idx))/0.5; end, X=z.*(ub-lb)+lb;
end
function s=levy(dim,beta)
sig=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta); s=randn(1,dim)*sig./(abs(randn(1,dim)).^(1/beta)+eps);
end
function y=bound(x,lb,ub), y=min(max(x,lb),ub); end
function T=fkineDH(dh,q)
T=eye(4); for i=1:size(dh,1), a=dh(i,1); al=dh(i,2); d=dh(i,3); th=q(i)+dh(i,4); A=[cos(th) -sin(th)*cos(al) sin(th)*sin(al) a*cos(th); sin(th) cos(th)*cos(al) -cos(th)*sin(al) a*sin(th); 0 sin(al) cos(al) d; 0 0 0 1]; T=T*A; end
end
