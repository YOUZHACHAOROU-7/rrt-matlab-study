%% main_ipso_353_timeopt_demo.m
% Paper-oriented demo:
% Improved PSO + 3-5-3 polynomial interpolation for manipulator time-optimal trajectory planning
% No Robotics Toolbox is required.
% Run this file directly in MATLAB.

clear; clc; close all; rng(2026);

%% 1. Waypoints in joint space, rad. Four points: q0 -> q1 -> q2 -> q3
% This follows the common 3-5-3 polynomial structure used in Chinese journal papers.
qWay = deg2rad([
    0   -45   35    0    25    0;
    25  -25   45   20    10   15;
    45  -10   20   35   -15   25;
    70   20  -10   45   -35   40]);

% Joint limits, velocity limits and acceleration limits. Adjust later to match real robot.
qLim = deg2rad(repmat([-170 170],6,1));
vMax = deg2rad([90 90 90 120 120 180]);
aMax = deg2rad([180 180 180 240 240 360]);

% Simplified 6-DOF DH parameters: [a alpha d theta_offset]
% This is a generic manipulator model only used for end-effector trajectory visualization.
dh = [
    0       pi/2   0.1625  0;
   -0.425   0      0       0;
   -0.3922  0      0       0;
    0       pi/2   0.1333  0;
    0      -pi/2   0.0997  0;
    0       0      0.0996  0];

%% 2. Optimizer settings
lb = [0.25 0.25 0.25];       % lower bound of t1,t2,t3
ub = [4.50 4.50 4.50];       % upper bound of t1,t2,t3
opts.nPop = 40;
opts.maxIter = 120;
opts.wMax = 0.90;
opts.wMin = 0.35;
opts.c1Max = 2.5;
opts.c1Min = 0.5;
opts.c2Max = 2.5;
opts.c2Min = 0.5;
opts.penalty = 1e4;
opts.levyBeta = 1.5;
opts.levyProb = 0.25;

obj = @(t) fitness_time_opt(t,qWay,qLim,vMax,aMax,opts.penalty);

%% 3. Run baseline PSO and improved Tent-Levy PSO
fprintf('Running baseline PSO...\n');
[psoBest,psoFit,psoCurve,psoInfo] = pso_basic(obj,lb,ub,opts);

fprintf('Running improved Tent-Levy PSO...\n');
[ipsoBest,ipsoFit,ipsoCurve,ipsoInfo] = pso_tent_levy(obj,lb,ub,opts);

%% 4. Evaluate optimized trajectories
psoTraj = build_353_trajectory(qWay,psoBest,250);
ipsoTraj = build_353_trajectory(qWay,ipsoBest,250);

fprintf('\n========== Result Summary ==========' );
fprintf('\nBaseline PSO best time: %.4f s, fitness: %.4f',sum(psoBest),psoFit);
fprintf('\nImproved PSO best time: %.4f s, fitness: %.4f',sum(ipsoBest),ipsoFit);
fprintf('\nTime reduction over baseline PSO: %.2f %%',100*(sum(psoBest)-sum(ipsoBest))/sum(psoBest));
fprintf('\nBaseline PSO t = [%.4f %.4f %.4f]',psoBest(1),psoBest(2),psoBest(3));
fprintf('\nImproved PSO t = [%.4f %.4f %.4f]\n',ipsoBest(1),ipsoBest(2),ipsoBest(3));

%% 5. Save table-style result
result.algorithm = {'PSO';'Tent-Levy-PSO'};
result.t1 = [psoBest(1);ipsoBest(1)];
result.t2 = [psoBest(2);ipsoBest(2)];
result.t3 = [psoBest(3);ipsoBest(3)];
result.totalTime = [sum(psoBest);sum(ipsoBest)];
result.bestFitness = [psoFit;ipsoFit];
T = struct2table(result);
disp(T);

%% 6. Plots for paper figures
figure('Color','w','Name','Convergence curve');
plot(psoCurve,'--','LineWidth',1.5); hold on;
plot(ipsoCurve,'-','LineWidth',2.0); grid on;
xlabel('Iteration'); ylabel('Best fitness');
legend('PSO','Tent-Levy-PSO','Location','northeast');
title('Convergence comparison');

figure('Color','w','Name','Joint displacement - improved PSO');
plot(ipsoTraj.time,rad2deg(ipsoTraj.q),'LineWidth',1.5); grid on;
xlabel('Time / s'); ylabel('Joint angle / deg');
title('Joint displacement curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');

figure('Color','w','Name','Joint velocity - improved PSO');
plot(ipsoTraj.time,rad2deg(ipsoTraj.qd),'LineWidth',1.5); grid on;
xlabel('Time / s'); ylabel('Joint velocity / deg/s');
title('Joint velocity curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');

figure('Color','w','Name','Joint acceleration - improved PSO');
plot(ipsoTraj.time,rad2deg(ipsoTraj.qdd),'LineWidth',1.5); grid on;
xlabel('Time / s'); ylabel('Joint acceleration / deg/s^2');
title('Joint acceleration curves');
legend('q1','q2','q3','q4','q5','q6','Location','best');

figure('Color','w','Name','End-effector path');
P = zeros(size(ipsoTraj.q,1),3);
for i = 1:size(ipsoTraj.q,1)
    T06 = fkine_dh(dh,ipsoTraj.q(i,:));
    P(i,:) = T06(1:3,4)';
end
plot3(P(:,1),P(:,2),P(:,3),'LineWidth',2); grid on; axis equal;
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m');
title('End-effector trajectory');

%% ================= local functions =================
function f = fitness_time_opt(t,qWay,qLim,vMax,aMax,penalty)
    if any(~isfinite(t)) || any(t<=0)
        f = 1e12; return;
    end
    traj = build_353_trajectory(qWay,t,180);
    totalTime = sum(t);
    vViol = max(0,abs(traj.qd) - vMax);
    aViol = max(0,abs(traj.qdd) - aMax);
    qLowViol = max(0,qLim(:,1)' - traj.q);
    qHighViol = max(0,traj.q - qLim(:,2)');
    smoothPenalty = sum(sum(abs(diff(traj.qdd,1,1)))) * 1e-4;
    f = totalTime + penalty*(sum(vViol(:).^2)+sum(aViol(:).^2)+sum(qLowViol(:).^2)+sum(qHighViol(:).^2)) + smoothPenalty;
end

function traj = build_353_trajectory(qWay,tSeg,nEach)
    % qWay: 4 x nJoint. tSeg: [t1,t2,t3].
    nJ = size(qWay,2);
    q = []; qd = []; qdd = []; time = [];
    tOffset = 0;
    for seg = 1:3
        if seg < 3
            tt = linspace(0,tSeg(seg),nEach)';
        else
            tt = linspace(0,tSeg(seg),nEach+1)';
        end
        if seg > 1
            tt = tt(2:end);
        end
        qSeg = zeros(numel(tt),nJ); qdSeg = qSeg; qddSeg = qSeg;
        for j = 1:nJ
            coeff = poly353_coeff_one_joint(qWay(:,j),tSeg);
            if seg == 1
                c = coeff(1:4);       % cubic: a13 a12 a11 a10
                qSeg(:,j) = polyval(c,tt);
                dc = polyder(c); ddc = polyder(dc);
            elseif seg == 2
                c = coeff(5:10);      % quintic: a25 ... a20
                qSeg(:,j) = polyval(c,tt);
                dc = polyder(c); ddc = polyder(dc);
            else
                c = coeff(11:14);     % cubic: a33 a32 a31 a30
                qSeg(:,j) = polyval(c,tt);
                dc = polyder(c); ddc = polyder(dc);
            end
            qdSeg(:,j) = polyval(dc,tt);
            qddSeg(:,j) = polyval(ddc,tt);
        end
        q = [q;qSeg]; qd = [qd;qdSeg]; qdd = [qdd;qddSeg]; %#ok<AGROW>
        time = [time; tOffset + tt]; %#ok<AGROW>
        tOffset = tOffset + tSeg(seg);
    end
    traj.q = q; traj.qd = qd; traj.qdd = qdd; traj.time = time;
end

function coeff = poly353_coeff_one_joint(theta,t)
    % Solve 14 coefficients for cubic-quintic-cubic interpolation.
    t1=t(1); t2=t(2); t3=t(3);
    th0=theta(1); th1=theta(2); th2=theta(3); th3=theta(4);
    A = zeros(14,14); b = zeros(14,1); row=0;
    % p1(0)=th0, p1'(0)=0, p1''(0)=0
    row=row+1; A(row,1:4) = poly_basis(3,0,0); b(row)=th0;
    row=row+1; A(row,1:4) = poly_basis(3,0,1); b(row)=0;
    row=row+1; A(row,1:4) = poly_basis(3,0,2); b(row)=0;
    % p1(t1)=th1, p2(0)=th1
    row=row+1; A(row,1:4) = poly_basis(3,t1,0); b(row)=th1;
    row=row+1; A(row,5:10)= poly_basis(5,0,0); b(row)=th1;
    % p1'(t1)=p2'(0), p1''(t1)=p2''(0)
    row=row+1; A(row,1:4)=poly_basis(3,t1,1); A(row,5:10)=-poly_basis(5,0,1);
    row=row+1; A(row,1:4)=poly_basis(3,t1,2); A(row,5:10)=-poly_basis(5,0,2);
    % p2(t2)=th2, p3(0)=th2
    row=row+1; A(row,5:10)=poly_basis(5,t2,0); b(row)=th2;
    row=row+1; A(row,11:14)=poly_basis(3,0,0); b(row)=th2;
    % p2'(t2)=p3'(0), p2''(t2)=p3''(0)
    row=row+1; A(row,5:10)=poly_basis(5,t2,1); A(row,11:14)=-poly_basis(3,0,1);
    row=row+1; A(row,5:10)=poly_basis(5,t2,2); A(row,11:14)=-poly_basis(3,0,2);
    % p3(t3)=th3, p3'(t3)=0, p3''(t3)=0
    row=row+1; A(row,11:14)=poly_basis(3,t3,0); b(row)=th3;
    row=row+1; A(row,11:14)=poly_basis(3,t3,1); b(row)=0;
    row=row+1; A(row,11:14)=poly_basis(3,t3,2); b(row)=0;
    coeff = (A\b)';
end

function v = poly_basis(order,t,der)
    powers = order:-1:0;
    v = zeros(1,order+1);
    for k = 1:numel(powers)
        p = powers(k);
        if p < der
            v(k) = 0;
        else
            mult = 1;
            for d = 0:der-1
                mult = mult*(p-d);
            end
            v(k) = mult*t^(p-der);
        end
    end
end

function [gbest,gfit,curve,info] = pso_basic(obj,lb,ub,opts)
    dim = numel(lb); nPop=opts.nPop; maxIter=opts.maxIter;
    X = rand(nPop,dim).*(ub-lb)+lb;
    V = zeros(nPop,dim);
    P = X; Pfit = arrayfun(@(i)obj(X(i,:)),1:nPop)';
    [gfit,idx] = min(Pfit); gbest = P(idx,:);
    curve = zeros(maxIter,1);
    for it=1:maxIter
        w = opts.wMax - (opts.wMax-opts.wMin)*it/maxIter;
        c1 = 2.0; c2 = 2.0;
        for i=1:nPop
            V(i,:) = w*V(i,:) + c1*rand(1,dim).*(P(i,:)-X(i,:)) + c2*rand(1,dim).*(gbest-X(i,:));
            X(i,:) = min(max(X(i,:) + V(i,:),lb),ub);
            fit = obj(X(i,:));
            if fit < Pfit(i)
                P(i,:) = X(i,:); Pfit(i)=fit;
            end
            if fit < gfit
                gbest = X(i,:); gfit=fit;
            end
        end
        curve(it)=gfit;
    end
    info.P = P; info.Pfit = Pfit;
end

function [gbest,gfit,curve,info] = pso_tent_levy(obj,lb,ub,opts)
    dim = numel(lb); nPop=opts.nPop; maxIter=opts.maxIter;
    X = tent_init(nPop,dim,lb,ub);
    V = zeros(nPop,dim);
    P = X; Pfit = arrayfun(@(i)obj(X(i,:)),1:nPop)';
    [gfit,idx] = min(Pfit); gbest = P(idx,:);
    curve = zeros(maxIter,1); stall = 0; lastBest = gfit;
    for it=1:maxIter
        w = opts.wMin + (opts.wMax-opts.wMin)*cos(pi*it/(2*maxIter))^2;
        c1 = opts.c1Max - (opts.c1Max-opts.c1Min)*(it/maxIter);
        c2 = opts.c2Min + (opts.c2Max-opts.c2Min)*(it/maxIter);
        for i=1:nPop
            V(i,:) = w*V(i,:) + c1*rand(1,dim).*(P(i,:)-X(i,:)) + c2*rand(1,dim).*(gbest-X(i,:));
            Xnew = X(i,:) + V(i,:);
            if rand < opts.levyProb
                step = levy_flight(dim,opts.levyBeta).*(X(i,:)-gbest);
                Xnew = Xnew + 0.08*step;
            end
            Xnew = min(max(Xnew,lb),ub);
            fit = obj(Xnew);
            % greedy selection
            if fit < Pfit(i)
                X(i,:) = Xnew; P(i,:) = Xnew; Pfit(i)=fit;
            else
                X(i,:) = Xnew;
            end
            if fit < gfit
                gbest = Xnew; gfit=fit;
            end
        end
        % small restart for stagnation
        if abs(lastBest-gfit) < 1e-8
            stall = stall + 1;
        else
            stall = 0; lastBest = gfit;
        end
        if stall > 15
            nRestart = max(2,round(0.15*nPop));
            [~,worstIdx] = maxk(Pfit,nRestart);
            X(worstIdx,:) = tent_init(nRestart,dim,lb,ub);
            for k=1:nRestart
                ii = worstIdx(k); P(ii,:) = X(ii,:); Pfit(ii)=obj(X(ii,:));
            end
            [tmpFit,tmpIdx] = min(Pfit);
            if tmpFit < gfit, gfit = tmpFit; gbest = P(tmpIdx,:); end
            stall = 0;
        end
        curve(it)=gfit;
    end
    info.P = P; info.Pfit = Pfit;
end

function X = tent_init(nPop,dim,lb,ub)
    z = rand(nPop,dim);
    for k = 1:8
        z = tent_map(z);
    end
    X = z.*(ub-lb)+lb;
end

function z = tent_map(z)
    mu = 0.499;
    idx = z < mu;
    z(idx) = z(idx)/mu;
    z(~idx) = (1-z(~idx))/(1-mu);
end

function step = levy_flight(dim,beta)
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1,dim)*sigma;
    v = randn(1,dim);
    step = u./(abs(v).^(1/beta)+eps);
end

function T = fkine_dh(dh,q)
    T = eye(4);
    for i=1:size(dh,1)
        a = dh(i,1); alpha = dh(i,2); d = dh(i,3); th = q(i)+dh(i,4);
        A = [cos(th), -sin(th)*cos(alpha),  sin(th)*sin(alpha), a*cos(th);
             sin(th),  cos(th)*cos(alpha), -cos(th)*sin(alpha), a*sin(th);
             0,        sin(alpha),          cos(alpha),         d;
             0,        0,                   0,                  1];
        T = T*A;
    end
end
