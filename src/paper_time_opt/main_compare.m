%% main_compare.m
% Paper-oriented demo: basic PSO vs improved PSO vs GWO for 3-5-3 polynomial
% time-optimal manipulator trajectory planning.
%
% Run this file directly in MATLAB.
% No external toolbox is required for the core optimization and trajectory plots.

clear; clc; close all;
rng(2026);

%% 1. Problem setup
prob = robot_dh_ur5e_problem();

% Four joint-space path points for 6-DOF manipulator, unit: rad.
% You can replace these points later with inverse-kinematic results from a real robot task.
prob.qPts = deg2rad([
     0   -55    60   -30    45    10;
    25   -40    45   -10    35   -15;
    50   -25    20    15    10   -30;
    70   -15     5    30   -10   -45]);

% Search variables: [t1, t2, t3].
prob.lb = [0.20 0.20 0.20];
prob.ub = [4.00 4.00 4.00];

% Joint limits for constraint penalty.
prob.qMin = deg2rad([-360 -360 -360 -360 -360 -360]);
prob.qMax = deg2rad([ 360  360  360  360  360  360]);
prob.vMax = deg2rad([120 120 150 180 180 240]);
prob.aMax = deg2rad([240 240 300 360 360 480]);
prob.nSamplesPerSeg = 80;
prob.penaltyWeight = 1e4;

% Optimizer parameters
opt.nPop = 40;
opt.maxIter = 120;
opt.dim = 3;
opt.lb = prob.lb;
opt.ub = prob.ub;
obj = @(x) fitness_time353(x, prob);

%% 2. Run algorithms
fprintf('Running basic PSO...\n');
resPSO = pso_basic(obj, opt);

fprintf('Running improved PSO...\n');
resIPSO = ipso_improved(obj, opt);

fprintf('Running basic GWO...\n');
resGWO = gwo_basic(obj, opt);

%% 3. Evaluate best trajectories
outPSO = evaluate353(resPSO.bestX, prob);
outIPSO = evaluate353(resIPSO.bestX, prob);
outGWO = evaluate353(resGWO.bestX, prob);

%% 4. Print table
fprintf('\n==================== Result Summary ====================\n');
fprintf('%-12s %-12s %-12s %-12s %-12s\n','Algorithm','t1','t2','t3','Fitness');
fprintf('%-12s %-12.4f %-12.4f %-12.4f %-12.4f\n','PSO',  resPSO.bestX(1),  resPSO.bestX(2),  resPSO.bestX(3),  resPSO.bestF);
fprintf('%-12s %-12.4f %-12.4f %-12.4f %-12.4f\n','IPSO', resIPSO.bestX(1), resIPSO.bestX(2), resIPSO.bestX(3), resIPSO.bestF);
fprintf('%-12s %-12.4f %-12.4f %-12.4f %-12.4f\n','GWO',  resGWO.bestX(1),  resGWO.bestX(2),  resGWO.bestX(3),  resGWO.bestF);
fprintf('========================================================\n');

%% 5. Plot convergence
figure('Color','w','Name','Convergence curves');
plot(resPSO.curve,'LineWidth',1.5); hold on;
plot(resIPSO.curve,'LineWidth',1.8);
plot(resGWO.curve,'LineWidth',1.5);
grid on; xlabel('Iteration'); ylabel('Best fitness');
title('Convergence comparison');
legend('PSO','Improved PSO','GWO','Location','northeast');

%% 6. Plot best joint trajectories from IPSO
figure('Color','w','Name','Joint position - IPSO');
plot(outIPSO.t, outIPSO.q, 'LineWidth',1.2); grid on;
xlabel('Time / s'); ylabel('Joint angle / rad');
title('Joint displacement curves optimized by improved PSO');
legend('q1','q2','q3','q4','q5','q6','Location','bestoutside');

figure('Color','w','Name','Joint velocity - IPSO');
plot(outIPSO.t, outIPSO.v, 'LineWidth',1.2); grid on;
xlabel('Time / s'); ylabel('Joint velocity / rad/s');
title('Joint velocity curves optimized by improved PSO');
legend('q1','q2','q3','q4','q5','q6','Location','bestoutside');

figure('Color','w','Name','Joint acceleration - IPSO');
plot(outIPSO.t, outIPSO.a, 'LineWidth',1.2); grid on;
xlabel('Time / s'); ylabel('Joint acceleration / rad/s^2');
title('Joint acceleration curves optimized by improved PSO');
legend('q1','q2','q3','q4','q5','q6','Location','bestoutside');

%% 7. Plot end-effector trajectory from IPSO
p = zeros(length(outIPSO.t),3);
for i = 1:length(outIPSO.t)
    T = fkine_dh(prob.dh, outIPSO.q(i,:));
    p(i,:) = T(1:3,4)';
end
figure('Color','w','Name','End-effector trajectory - IPSO');
plot3(p(:,1), p(:,2), p(:,3), 'LineWidth',1.8); grid on; axis equal;
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m');
title('End-effector trajectory optimized by improved PSO');

%% 8. Save results for paper tables
if ~exist('results','dir'), mkdir('results'); end
T = table({'PSO';'IPSO';'GWO'}, ...
    [resPSO.bestX(1);resIPSO.bestX(1);resGWO.bestX(1)], ...
    [resPSO.bestX(2);resIPSO.bestX(2);resGWO.bestX(2)], ...
    [resPSO.bestX(3);resIPSO.bestX(3);resGWO.bestX(3)], ...
    [sum(resPSO.bestX);sum(resIPSO.bestX);sum(resGWO.bestX)], ...
    [resPSO.bestF;resIPSO.bestF;resGWO.bestF], ...
    'VariableNames',{'Algorithm','t1','t2','t3','TotalTime','Fitness'});
writetable(T, fullfile('results','comparison_results.csv'));
save(fullfile('results','paper_demo_results.mat'),'resPSO','resIPSO','resGWO','outPSO','outIPSO','outGWO','prob');
fprintf('\nSaved: results/comparison_results.csv and results/paper_demo_results.mat\n');
