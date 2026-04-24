%% main_paper_demo.m
% 单文件演示：改进RRT + 路径剪枝 + 时间参数化
clear;clc;close all; rng(7);

%% 1) 参数
jointLimits = repmat([-pi,pi],6,1);
qStart = deg2rad([-70 -30 35 10 20 0]);
qGoal  = deg2rad([45 25 -20 35 -30 15]);
dh = [0 pi/2 0.34 0; 0.30 0 0 0; 0.22 0 0 0; 0 pi/2 0.12 0; 0 -pi/2 0.08 0; 0 0 0.10 0];
obs(1).min = [0.10 -0.45 0.10]; obs(1).max = [0.35 -0.10 0.55];
obs(2).min = [0.15  0.05 0.25]; obs(2).max = [0.42  0.28 0.68];
obs(3).min = [-0.25 -0.08 0.18]; obs(3).max = [-0.05 0.18 0.62];
P.maxIter = 3500; P.stepMin = 0.10; P.stepMax = 0.40; P.goalBiasMin = 0.10; P.goalBiasMax = 0.45;
P.goalThreshold = 0.22; P.rewireRadius = 0.50; P.validationStep = 0.08; P.shortcutIters = 180; P.linkSampleCount = 12;
vMax = deg2rad([80 80 90 120 120 150]); aMax = deg2rad([120 120 140 180 180 220]);
assert(~stateCollision(qStart,dh,obs,P.linkSampleCount),'qStart collision');
assert(~stateCollision(qGoal,dh,obs,P.linkSampleCount),'qGoal collision');

%% 2) 改进RRT
fprintf('开始规划...\n'); tic;
[pathRaw,tree,stats] = improvedRRT(qStart,qGoal,jointLimits,dh,obs,P);
planT = toc;
if isempty(pathRaw), error('未找到可行路径'); end
pathOpt = shortcutPath(pathRaw,dh,obs,P);
traj = timeParameterize(pathOpt,vMax,aMax);

%% 3) 指标
fprintf('\n===== 结果摘要 =====\n');
fprintf('原始路径长度: %.4f rad\n', pathLength(pathRaw));
fprintf('平滑路径长度: %.4f rad\n', pathLength(pathOpt));
fprintf('原始节点数  : %d\n', size(pathRaw,1));
fprintf('平滑节点数  : %d\n', size(pathOpt,1));
fprintf('规划耗时    : %.4f s\n', planT);
fprintf('执行总时间  : %.4f s\n', traj.totalTime);
fprintf('扩展成功次数: %d\n', stats.expandCount);

%% 4) 可视化
figure('Color','w','Name','Workspace');
subplot(1,2,1); showScene(dh,obs,qStart,qGoal,pathRaw,'原始路径');
subplot(1,2,2); showScene(dh,obs,qStart,qGoal,pathOpt,'平滑后路径');
figure('Color','w','Name','Joint Path');
plot(pathRaw,'--','LineWidth',1.0); hold on; plot(pathOpt,'LineWidth',1.8);
grid on; xlabel('路径点'); ylabel('关节角/rad'); title('关节空间路径对比');
figure('Color','w','Name','Time Allocation');
stairs(traj.tBreaks,[traj.segTime;traj.segTime(end)],'LineWidth',1.8); grid on;
xlabel('累计时间/s'); ylabel('分段时间/s'); title('时间参数化结果');

%% ====== functions ======
function [path,tree,stats] = improvedRRT(qStart,qGoal,jointLimits,dh,obs,P)
tree.nodes = qStart; tree.parent = 0; tree.cost = 0; initDist = norm(qGoal-qStart); bestDist = initDist; goalNode = [];
stats.expandCount = 0; stats.bestDistance = initDist;
for iter = 1:P.maxIter
    progress = 1 - bestDist/max(initDist,eps);
    pGoal = P.goalBiasMin + progress*(P.goalBiasMax-P.goalBiasMin);
    if rand < pGoal, qRand = qGoal; else, qRand = jointLimits(:,1)' + rand(1,6).*(jointLimits(:,2)'-jointLimits(:,1)'); end
    idxNear = nearestNode(tree.nodes,qRand); qNear = tree.nodes(idxNear,:);
    step = P.stepMin + min(norm(qGoal-qNear)/max(initDist,eps),1)*(P.stepMax-P.stepMin);
    qNew = steer(qNear,qRand,step);
    if ~transitionValid(qNear,qNew,dh,obs,P), continue; end
    nearIdx = findNearNodes(tree.nodes,qNew,P.rewireRadius);
    [parentIdx,parentCost] = chooseParent(tree,nearIdx,idxNear,qNew,dh,obs,P);
    tree.nodes(end+1,:) = qNew; tree.parent(end+1,1) = parentIdx; tree.cost(end+1,1) = parentCost; newIdx = size(tree.nodes,1);
    stats.expandCount = stats.expandCount + 1; tree = rewire(tree,nearIdx,newIdx,dh,obs,P);
    dGoal = norm(qNew-qGoal); if dGoal < bestDist, bestDist = dGoal; stats.bestDistance = dGoal; end
    if dGoal < P.goalThreshold && transitionValid(qNew,qGoal,dh,obs,P)
        tree.nodes(end+1,:) = qGoal; tree.parent(end+1,1) = newIdx; tree.cost(end+1,1) = tree.cost(newIdx)+norm(qGoal-qNew); goalNode = size(tree.nodes,1); break;
    end
end
if isempty(goalNode), path = []; else, path = backtrack(tree,goalNode); end
end
function idx = nearestNode(nodes,q), [~,idx] = min(sum((nodes-q).^2,2)); end
function qNew = steer(q1,q2,step), v = q2-q1; d = norm(v); if d<1e-12, qNew=q1; elseif d<=step, qNew=q2; else, qNew=q1+step*v/d; end, end
function idxs = findNearNodes(nodes,q,r), idxs = find(sqrt(sum((nodes-q).^2,2))<=r); end
function [bestIdx,bestCost] = chooseParent(tree,nearIdx,fallbackIdx,qNew,dh,obs,P)
bestIdx = fallbackIdx; bestCost = tree.cost(fallbackIdx)+norm(tree.nodes(fallbackIdx,:)-qNew);
for k=1:numel(nearIdx)
    idx = nearIdx(k); c = tree.cost(idx)+norm(tree.nodes(idx,:)-qNew);
    if c<bestCost && transitionValid(tree.nodes(idx,:),qNew,dh,obs,P), bestIdx = idx; bestCost = c; end
end
end
function tree = rewire(tree,nearIdx,newIdx,dh,obs,P)
qNew = tree.nodes(newIdx,:); cNew = tree.cost(newIdx);
for k=1:numel(nearIdx)
    idx = nearIdx(k); if idx==newIdx || idx==tree.parent(newIdx), continue; end
    cTry = cNew + norm(tree.nodes(idx,:)-qNew);
    if cTry<tree.cost(idx) && transitionValid(qNew,tree.nodes(idx,:),dh,obs,P), tree.parent(idx)=newIdx; tree.cost(idx)=cTry; end
end
end
function ok = transitionValid(qA,qB,dh,obs,P)
n = max(2,ceil(norm(qB-qA)/P.validationStep)); ok = true;
for i=0:n
    q = qA + (i/n)*(qB-qA);
    if stateCollision(q,dh,obs,P.linkSampleCount), ok = false; return; end
end
end
function flag = stateCollision(q,dh,obs,sampleCount)
pts = jointPoints(q,dh); flag = false;
for i=1:size(pts,1)-1
    p1 = pts(i,:); p2 = pts(i+1,:);
    for s = linspace(0,1,sampleCount)
        p = p1 + s*(p2-p1);
        for j=1:numel(obs)
            if all(p>=obs(j).min) && all(p<=obs(j).max), flag = true; return; end
        end
    end
end
end
function pts = jointPoints(q,dh)
T = eye(4); pts = zeros(size(dh,1)+1,3); pts(1,:) = [0 0 0];
for i=1:size(dh,1)
    a=dh(i,1); alpha=dh(i,2); d=dh(i,3); theta=q(i)+dh(i,4);
    A=[cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta); sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); 0 sin(alpha) cos(alpha) d; 0 0 0 1];
    T=T*A; pts(i+1,:)=T(1:3,4)';
end
end
function path = backtrack(tree,idx)
path = tree.nodes(idx,:); while tree.parent(idx)~=0, idx = tree.parent(idx); path = [tree.nodes(idx,:); path]; end
end
function newPath = shortcutPath(path,dh,obs,P)
newPath = path; if size(path,1)<=2, return; end
for it = 1:P.shortcutIters
    if size(newPath,1)<=2, break; end
    i = randi([1,size(newPath,1)-1]); j = randi([i+1,size(newPath,1)]); if j<=i+1, continue; end
    if transitionValid(newPath(i,:),newPath(j,:),dh,obs,P), newPath = [newPath(1:i,:); newPath(j:end,:)]; end
end
end
function L = pathLength(path), if size(path,1)<2, L=0; else, dq = diff(path,1,1); L = sum(sqrt(sum(dq.^2,2))); end, end
function traj = timeParameterize(path,vMax,aMax)
nSeg = size(path,1)-1; segTime = zeros(nSeg,1);
for i=1:nSeg
    dq = abs(path(i+1,:)-path(i,:)); tJ = zeros(1,numel(dq));
    for j=1:numel(dq)
        tAcc = vMax(j)/aMax(j); dAcc = 0.5*aMax(j)*tAcc^2;
        if dq(j)<=2*dAcc, tJ(j)=2*sqrt(dq(j)/aMax(j)); else, tJ(j)=2*tAcc + (dq(j)-2*dAcc)/vMax(j); end
    end
    segTime(i)=max(tJ);
end
traj.segTime = segTime; traj.tBreaks = [0;cumsum(segTime)]; traj.totalTime = sum(segTime);
end
function showScene(dh,obs,qStart,qGoal,path,ttl)
hold on; axis equal; grid on; view(3); xlabel('X/m'); ylabel('Y/m'); zlabel('Z/m'); title(ttl);
for i=1:numel(obs), drawBox(obs(i).min,obs(i).max,[0.85 0.3 0.3],0.25); end
plotRobot(jointPoints(qStart,dh),[0.1 0.4 0.9],2.4); plotRobot(jointPoints(qGoal,dh),[0.1 0.7 0.2],2.4);
pickIdx = unique(round(linspace(1,size(path,1),min(12,size(path,1)))));
for k=1:numel(pickIdx), plotRobot(jointPoints(path(pickIdx(k),:),dh),[0.35 0.35 0.35],1.0); end
end
function plotRobot(pts,c,lw), plot3(pts(:,1),pts(:,2),pts(:,3),'-o','Color',c,'LineWidth',lw,'MarkerSize',4,'MarkerFaceColor',c); end
function drawBox(bmin,bmax,c,a)
v=[bmin(1) bmin(2) bmin(3); bmax(1) bmin(2) bmin(3); bmax(1) bmax(2) bmin(3); bmin(1) bmax(2) bmin(3); bmin(1) bmin(2) bmax(3); bmax(1) bmin(2) bmax(3); bmax(1) bmax(2) bmax(3); bmin(1) bmax(2) bmax(3)];
f=[1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8]; patch('Vertices',v,'Faces',f,'FaceColor',c,'FaceAlpha',a,'EdgeColor',c,'LineWidth',1.0);
end
