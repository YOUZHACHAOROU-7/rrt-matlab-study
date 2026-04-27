%% main_paper_18d_full_experiment.m
% Paper-level 18-dimensional experiment.
% Each joint owns [t1,t2,t3]. Objective is the average motion time plus constraints.
% Algorithms: Fixed-353, PSO, GWO, WOA, Tent-Levy-PSO.
% Run from repository root: run('src/main_paper_18d_full_experiment.m')

clear; clc; close all; rng(2026);

rootDir = pwd;
if endsWith(rootDir, [filesep 'src'])
    rootDir = fileparts(rootDir);
end
outDir = fullfile(rootDir,'results_18d_full');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% 1. Waypoints: 4 points x 6 joints, rad
qWay = deg2rad([
     0  -45   35    0   25    0;
    25  -25   45   20   10   15;
    45  -10   20   35  -15   25;
    70   20  -10   45  -35   40]);

qLim = deg2rad(repmat([-170 170],6,1));
vMax = deg2rad([90 90 90 120 120 180]);
aMax = deg2rad([180 180 180 240 240 360]);
dh = [0 pi/2 0.1625 0; -0.425 0 0 0; -0.3922 0 0 0; 0 pi/2 0.1333 0; 0 -pi/2 0.0997 0; 0 0 0.0996 0];

%% 2. Settings
nJoint = 6; dim = 18;
lb = repmat([0.25 0.25 0.25],1,nJoint);
ub = repmat([4.50 4.50 4.50],1,nJoint);
fixedX = repmat([2.5 2.5 2.5],1,nJoint);

opt.nPop = 45;
opt.maxIter = 180;
opt.wMax = 0.90; opt.wMin = 0.35;
opt.c1Max = 2.5; opt.c1Min = 0.5;
opt.c2Max = 2.5; opt.c2Min = 0.5;
opt.penalty = 1e4;
opt.levyProb = 0.35;
opt.levyBeta = 1.5;
opt.dim = dim;

obj = @(x) fitness18D(x,qWay,qLim,vMax,aMax,opt.penalty);
fixedF = obj(fixedX); fixedTime = mean(reshape(fixedX,3,6),'all')*3;

algNames = {'PSO','GWO','WOA','Tent-Levy-PSO'};
nAlg = numel(algNames);
nRun = 20;
bestX = zeros(nRun,dim,nAlg);
bestF = zeros(nRun,nAlg);
bestTime = zeros(nRun,nAlg);
curves = zeros(opt.maxIter,nAlg,nRun);

fprintf('18D full experiment starts: %d repeated runs.\n', nRun);
for r = 1:nRun
    rng(9000+r);
    [x1,f1,c1] = basicPSO(obj,lb,ub,opt);
    [x2,f2,c2] = gwoBasic(obj,lb,ub,opt);
    [x3,f3,c3] = woaBasic(obj,lb,ub,opt);
    [x4,f4,c4] = tlpso(obj,lb,ub,opt);
    Xs = cat(3,x1,x2,x3,x4); Fs=[f1 f2 f3 f4]; Cs={c1,c2,c3,c4};
    for a=1:nAlg
        bestX(r,:,a)=Xs(:,:,a);
        bestF(r,a)=Fs(a);
        bestTime(r,a)=meanJointTotalTime(bestX(r,:,a));
        curves(:,a,r)=Cs{a};
    end
    fprintf('Run %02d/%02d | PSO %.4f | GWO %.4f | WOA %.4f | TLPSO %.4f\n', ...
        r,nRun,bestTime(r,1),bestTime(r,2),bestTime(r,3),bestTime(r,4));
end

%% 3. Summary
FixedRow = table({'Fixed-353'},fixedTime,fixedF,NaN,NaN,NaN, ...
    'VariableNames',{'Algorithm','BestTime','BestFitness','MeanTime','StdTime','MeanFitness'});
Rows = FixedRow;
for a=1:nAlg
    [bf,idx] = min(bestF(:,a));
    bt = bestTime(idx,a);
    row = table(algNames(a),bt,bf,mean(bestTime(:,a)),std(bestTime(:,a)),mean(bestF(:,a)), ...
        'VariableNames',{'Algorithm','BestTime','BestFitness','MeanTime','StdTime','MeanFitness'});
    Rows = [Rows; row]; %#ok<AGROW>
end
Summary = Rows;
disp('========== 18D Paper Summary Table ==========');
disp(Summary);
writetable(Summary,fullfile(outDir,'paper_summary_18d.csv'));

ReductionVsFixed = 100*(fixedTime-Summary.BestTime)/fixedTime;
ReductionTable = table(Summary.Algorithm,Summary.BestTime,ReductionVsFixed, ...
    'VariableNames',{'Algorithm','BestTime','ReductionComparedWithFixedPercent'});
disp('========== Reduction Compared With Fixed 3-5-3 ==========');
disp(ReductionTable);
writetable(ReductionTable,fullfile(outDir,'reduction_vs_fixed_18d.csv'));

Detail = table((1:nRun)',bestTime(:,1),bestF(:,1),bestTime(:,2),bestF(:,2),bestTime(:,3),bestF(:,3),bestTime(:,4),bestF(:,4), ...
    'VariableNames',{'RunID','PSO_Time','PSO_Fit','GWO_Time','GWO_Fit','WOA_Time','WOA_Fit','TLPSO_Time','TLPSO_Fit'});
writetable(Detail,fullfile(outDir,'repeated_runs_18d.csv'));

%% 4. Plots
fig1=figure('Color','w','Name','18D average convergence'); hold on;
styles={'--',':','-.','-'};
for a=1:nAlg
    plot(mean(curves(:,a,:),3),styles{a},'LineWidth',1.8);
end
grid on; xlabel('Iteration'); ylabel('Mean best fitness');
legend(algNames,'Location','northeast'); title('Average convergence curves');
saveas(fig1,fullfile(outDir,'fig_18d_average_convergence.png'));

[~,bestRun] = min(bestF(:,4));
bestTL = bestX(bestRun,:,4);
traj = buildTraj18D(qWay,bestTL,260);

fig2=figure('Color','w','Name','18D joint displacement');
plot(traj.time,rad2deg(traj.q),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint angle / deg'); title('Joint displacement curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig2,fullfile(outDir,'fig_18d_joint_displacement.png'));

fig3=figure('Color','w','Name','18D joint velocity');
plot(traj.time,rad2deg(traj.qd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint velocity / deg/s'); title('Joint velocity curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig3,fullfile(outDir,'fig_18d_joint_velocity.png'));

fig4=figure('Color','w','Name','18D joint acceleration');
plot(traj.time,rad2deg(traj.qdd),'LineWidth',1.4); grid on;
xlabel('Time / s'); ylabel('Joint acceleration / deg/s^2'); title('Joint acceleration curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');
saveas(fig4,fullfile(outDir,'fig_18d_joint_acceleration.png'));

P=zeros(size(traj.q,1),3);
for i=1:size(traj.q,1)
    T=fkineDH(dh,traj.q(i,:)); P(i,:)=T(1:3,4)';
end
fig5=figure('Color','w','Name','18D end-effector trajectory');
plot3(P(:,1),P(:,2),P(:,3),'LineWidth',2); grid on; axis equal;
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m'); title('End-effector trajectory');
saveas(fig5,fullfile(outDir,'fig_18d_end_effector_trajectory.png'));

fprintf('\nAll 18D tables and figures have been saved to: %s\n',outDir);

%% ================= local functions =================
function f = fitness18D(x,qWay,qLim,vMax,aMax,penalty)
    if any(x<=0)||any(~isfinite(x)), f=1e12; return; end
    tMat=reshape(x,3,6); tr=buildTraj18D(qWay,x,150);
    motionTime = mean(sum(tMat,1));
    syncPenalty = 0.04*std(sum(tMat,1));
    vViol=max(0,abs(tr.qd)-vMax);
    aViol=max(0,abs(tr.qdd)-aMax);
    qViol=max(0,qLim(:,1)'-tr.q)+max(0,tr.q-qLim(:,2)');
    smooth=sum(abs(diff(tr.qdd,1,1)),'all')*1e-4;
    f=motionTime+syncPenalty+penalty*(sum(vViol.^2,'all')+sum(aViol.^2,'all')+sum(qViol.^2,'all'))+smooth;
end

function T=meanJointTotalTime(x)
    T=mean(sum(reshape(x,3,6),1));
end

function tr=buildTraj18D(qWay,x,nEach)
    tMat=reshape(x,3,6); totalT=mean(sum(tMat,1)); N=3*nEach-2;
    q=zeros(N,6); qd=q; qdd=q; time=linspace(0,totalT,N)';
    for j=1:6
        trj=buildOneJoint(qWay(:,j),tMat(:,j)',nEach);
        tq=linspace(0,sum(tMat(:,j)),numel(trj.t))';
        q(:,j)=interp1(tq,trj.q,time,'pchip','extrap');
        qd(:,j)=interp1(tq,trj.qd,time,'pchip','extrap');
        qdd(:,j)=interp1(tq,trj.qdd,time,'pchip','extrap');
    end
    tr.q=q; tr.qd=qd; tr.qdd=qdd; tr.time=time;
end

function tr=buildOneJoint(theta,tSeg,nEach)
    q=[]; qd=[]; qdd=[]; time=[]; offset=0; coef=coef353(theta,tSeg);
    for s=1:3
        tt=linspace(0,tSeg(s),nEach)'; if s>1, tt=tt(2:end); end
        if s==1, c=coef(1:4); elseif s==2, c=coef(5:10); else, c=coef(11:14); end
        dc=polyder(c); ddc=polyder(dc);
        q=[q;polyval(c,tt)]; qd=[qd;polyval(dc,tt)]; qdd=[qdd;polyval(ddc,tt)]; time=[time;offset+tt]; offset=offset+tSeg(s);
    end
    tr.q=q; tr.qd=qd; tr.qdd=qdd; tr.t=time;
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

function [gb,gf,curve]=tlpso(obj,lb,ub,opt)
    dim=numel(lb); n=opt.nPop; X=tentInit(n,dim,lb,ub); V=zeros(n,dim); P=X; Pf=zeros(n,1);
    for i=1:n, Pf(i)=obj(X(i,:)); end
    [gf,id]=min(Pf); gb=P(id,:); curve=zeros(opt.maxIter,1); stall=0; last=gf;
    for it=1:opt.maxIter
        w=opt.wMin+(opt.wMax-opt.wMin)*cos(pi*it/(2*opt.maxIter))^2;
        c1=opt.c1Max-(opt.c1Max-opt.c1Min)*it/opt.maxIter;
        c2=opt.c2Min+(opt.c2Max-opt.c2Min)*it/opt.maxIter;
        for i=1:n
            V(i,:)=w*V(i,:)+c1*rand(1,dim).*(P(i,:)-X(i,:))+c2*rand(1,dim).*(gb-X(i,:));
            Xn=X(i,:)+V(i,:);
            if rand<opt.levyProb, Xn=Xn+0.06*levy(dim,opt.levyBeta).*(X(i,:)-gb); end
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

function [Alpha_pos,Alpha_score,curve]=gwoBasic(obj,lb,ub,opt)
    dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb;
    Alpha_score=inf; Beta_score=inf; Delta_score=inf; Alpha_pos=zeros(1,dim); Beta_pos=Alpha_pos; Delta_pos=Alpha_pos;
    curve=zeros(opt.maxIter,1);
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
        for i=1:n
            r1=rand(1,dim); r2=rand(1,dim); A1=2*a*r1-a; C1=2*r2; D_alpha=abs(C1.*Alpha_pos-X(i,:)); X1=Alpha_pos-A1.*D_alpha;
            r1=rand(1,dim); r2=rand(1,dim); A2=2*a*r1-a; C2=2*r2; D_beta=abs(C2.*Beta_pos-X(i,:)); X2=Beta_pos-A2.*D_beta;
            r1=rand(1,dim); r2=rand(1,dim); A3=2*a*r1-a; C3=2*r2; D_delta=abs(C3.*Delta_pos-X(i,:)); X3=Delta_pos-A3.*D_delta;
            X(i,:)=(X1+X2+X3)/3;
        end
        curve(it)=Alpha_score;
    end
end

function [bestX,bestF,curve]=woaBasic(obj,lb,ub,opt)
    dim=numel(lb); n=opt.nPop; X=rand(n,dim).*(ub-lb)+lb; bestX=zeros(1,dim); bestF=inf; curve=zeros(opt.maxIter,1); b=1;
    for it=1:opt.maxIter
        for i=1:n
            X(i,:)=min(max(X(i,:),lb),ub); fit=obj(X(i,:));
            if fit<bestF, bestF=fit; bestX=X(i,:); end
        end
        a=2-2*it/opt.maxIter;
        for i=1:n
            r=rand; A=2*a*r-a; C=2*rand; p=rand; l=-1+2*rand;
            if p<0.5
                if abs(A)<1
                    D=abs(C*bestX-X(i,:)); X(i,:)=bestX-A.*D;
                else
                    xr=X(randi(n),:); D=abs(C*xr-X(i,:)); X(i,:)=xr-A.*D;
                end
            else
                D=abs(bestX-X(i,:)); X(i,:)=D.*exp(b*l).*cos(2*pi*l)+bestX;
            end
        end
        curve(it)=bestF;
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
