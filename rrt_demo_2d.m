clc;
clear;
close all;

%% 地图范围
x_max = 100;
y_max = 100;

%% 起点和终点
start = [10, 10];
goal  = [90, 90];

%% 障碍物 [x, y, w, h]
obstacles = [
    20, 20, 15, 40;
    50, 50, 20, 20;
    70, 10, 10, 50
];

%% 参数
max_iter = 2000;      % 最大迭代次数
step_size = 5;        % 每次扩展步长
goal_threshold = 8;   % 到目标点的判定距离

%% 树结构
nodes = start;        % 每一行是一个节点
parents = 0;          % 每个节点的父节点索引，第一个节点没有父节点

goal_reached = false;
goal_index = -1;

figure;
hold on;
axis([0 x_max 0 y_max]);
axis equal;
grid on;
title('2D RRT Demo');

% 画障碍物
for i = 1:size(obstacles,1)
    rectangle('Position', obstacles(i,:), 'FaceColor', [0.5 0.5 0.5]);
end

% 画起点终点
plot(start(1), start(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(goal(1), goal(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

for iter = 1:max_iter
    %% 1. 随机采样
    rand_point = [rand*x_max, rand*y_max];

    %% 2. 找最近节点
    dists = sqrt(sum((nodes - rand_point).^2, 2));
    [~, nearest_idx] = min(dists);
    nearest_node = nodes(nearest_idx, :);

    %% 3. 朝随机点扩展一步
    direction = rand_point - nearest_node;
    direction = direction / norm(direction);
    new_node = nearest_node + step_size * direction;

    %% 边界检查
    if new_node(1) < 0 || new_node(1) > x_max || new_node(2) < 0 || new_node(2) > y_max
        continue;
    end

    %% 4. 碰撞检测
    if is_collision(new_node, obstacles)
        continue;
    end

    %% 5. 加入树
    nodes = [nodes; new_node];
    parents = [parents; nearest_idx];

    plot([nearest_node(1), new_node(1)], [nearest_node(2), new_node(2)], 'b-');
    drawnow limitrate;

    %% 6. 判断是否到达目标
    if norm(new_node - goal) < goal_threshold
        goal_reached = true;
        goal_index = size(nodes,1);
        break;
    end
end

%% 7. 回溯路径
if goal_reached
    path = goal;
    current_idx = goal_index;

    while current_idx ~= 1
        path = [nodes(current_idx,:); path];
        current_idx = parents(current_idx);
    end
    path = [start; path];

    % 画最终路径
    plot(path(:,1), path(:,2), 'r-', 'LineWidth', 2);
    disp('找到路径！');
else
    disp('未找到路径。');
end

%% ===== 本地函数 =====
function flag = is_collision(point, obstacles)
    flag = false;
    for k = 1:size(obstacles,1)
        obs = obstacles(k,:);
        if point(1) >= obs(1) && point(1) <= obs(1)+obs(3) && ...
           point(2) >= obs(2) && point(2) <= obs(2)+obs(4)
            flag = true;
            return;
        end
    end
end