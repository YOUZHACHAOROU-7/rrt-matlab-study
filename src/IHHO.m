function [Rabbit_Energy, Rabbit_Location, Convergence_curve] = IHHO(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% IHHO  Multi-strategy Improved Harris Hawks Optimizer.
%
% Strategies:
% 1) Tent chaotic opposition-based initialization.
% 2) Nonlinear escape energy factor.
% 3) Stagnation-aware Levy-Cauchy mutation with greedy selection.
%
% Usage:
%   [bestScore,bestPos,curve] = IHHO(N,MaxIt,lb,ub,dim,fobj)

    if numel(lb)==1, lb = lb * ones(1,dim); end
    if numel(ub)==1, ub = ub * ones(1,dim); end

    % Strategy 1: Tent chaotic initialization + opposition-based selection.
    X0 = tentInit(SearchAgents_no, dim, lb, ub);
    Xopp = lb + ub - X0;
    Xall = [X0; Xopp];
    fitAll = inf(size(Xall,1),1);
    for i = 1:size(Xall,1)
        Xall(i,:) = bound(Xall(i,:),lb,ub);
        fitAll(i) = fobj(Xall(i,:));
    end
    [~, idx] = sort(fitAll,'ascend');
    X = Xall(idx(1:SearchAgents_no),:);

    Rabbit_Location = zeros(1,dim);
    Rabbit_Energy = inf;
    Convergence_curve = zeros(Max_iter,1);
    fit = inf(SearchAgents_no,1);
    stale = zeros(SearchAgents_no,1);

    for i = 1:SearchAgents_no
        fit(i) = fobj(X(i,:));
        if fit(i) < Rabbit_Energy
            Rabbit_Energy = fit(i);
            Rabbit_Location = X(i,:);
        end
    end

    eliteCount = max(1,round(0.08*SearchAgents_no));

    for t = 1:Max_iter
        tau = t/Max_iter;

        % Update rabbit and keep sorted elites.
        for i = 1:SearchAgents_no
            X(i,:) = bound(X(i,:),lb,ub);
            fit(i) = fobj(X(i,:));
            if fit(i) < Rabbit_Energy
                Rabbit_Energy = fit(i);
                Rabbit_Location = X(i,:);
            end
        end
        [fit, order] = sort(fit,'ascend');
        X = X(order,:);
        elites = X(1:eliteCount,:);
        eliteFit = fit(1:eliteCount);
        Rabbit_Energy = fit(1);
        Rabbit_Location = X(1,:);
        X_mean = mean(X,1);

        % Strategy 2: nonlinear escaping energy, smoother than original HHO.
        E1 = 2 * (1 - tau^1.7) * cos(pi*tau/2);

        X_old = X;
        fit_old = fit;

        for i = 1:SearchAgents_no
            E0 = 2*rand - 1;
            Escaping_Energy = E1 * E0;
            q = rand;
            r = rand;
            J = 2*(1-rand);

            if abs(Escaping_Energy) >= 1
                % Exploration phase.
                rand_idx = randi(SearchAgents_no);
                X_rand = X(rand_idx,:);
                if q < 0.5
                    X_new = X_rand - rand(1,dim).*abs(X_rand - 2*rand(1,dim).*X(i,:));
                else
                    X_new = (Rabbit_Location - X_mean) - rand(1,dim).*((ub-lb).*rand(1,dim) + lb);
                end
            else
                % Exploitation phase.
                if r >= 0.5 && abs(Escaping_Energy) < 0.5
                    X_new = Rabbit_Location - Escaping_Energy .* abs(Rabbit_Location - X(i,:));
                elseif r >= 0.5 && abs(Escaping_Energy) >= 0.5
                    X_new = (Rabbit_Location - X(i,:)) - Escaping_Energy .* abs(J*Rabbit_Location - X(i,:));
                elseif r < 0.5 && abs(Escaping_Energy) >= 0.5
                    Y = Rabbit_Location - Escaping_Energy .* abs(J*Rabbit_Location - X(i,:));
                    Y = bound(Y,lb,ub);
                    Z = Y + randn(1,dim).*levyStep(dim,1.5);
                    Z = bound(Z,lb,ub);
                    if fobj(Y) < fit(i)
                        X_new = Y;
                    elseif fobj(Z) < fit(i)
                        X_new = Z;
                    else
                        X_new = X(i,:);
                    end
                else
                    Y = Rabbit_Location - Escaping_Energy .* abs(J*Rabbit_Location - X_mean);
                    Y = bound(Y,lb,ub);
                    Z = Y + randn(1,dim).*levyStep(dim,1.5);
                    Z = bound(Z,lb,ub);
                    if fobj(Y) < fit(i)
                        X_new = Y;
                    elseif fobj(Z) < fit(i)
                        X_new = Z;
                    else
                        X_new = X(i,:);
                    end
                end
            end

            X_new = bound(X_new,lb,ub);
            f_new = fobj(X_new);

            % Strategy 3: trigger mutation mainly on stagnant individuals.
            isStagnant = stale(i) >= 4;
            pm = 0.18*(1-tau) + 0.03;
            if isStagnant || rand < pm
                levyScale = 0.025*(1-tau+0.05);
                cauchyScale = 0.008*(1-tau+0.05);  % lower late Cauchy disturbance for stability.
                localScale = 0.010*(1-tau+0.05);
                X_mut = X_new ...
                    + levyScale * levyStep(dim,1.5) .* (X_new - Rabbit_Location) ...
                    + cauchyScale * cauchyRand(1,dim) .* (ub-lb) ...
                    + localScale * randn(1,dim) .* abs(Rabbit_Location - X_new);
                X_mut = bound(X_mut,lb,ub);
                f_mut = fobj(X_mut);
                if f_mut < f_new
                    X_new = X_mut;
                    f_new = f_mut;
                end
            end

            % Greedy selection.
            if f_new < fit_old(i)
                X(i,:) = X_new;
                fit(i) = f_new;
                stale(i) = 0;
            else
                X(i,:) = X_old(i,:);
                fit(i) = fit_old(i);
                stale(i) = stale(i) + 1;
            end
        end

        % Elite preservation avoids losing the best solutions after disturbance.
        [~, worstIdx] = sort(fit,'descend');
        replaceIdx = worstIdx(1:eliteCount);
        X(replaceIdx,:) = elites;
        fit(replaceIdx) = eliteFit;

        [Rabbit_Energy, bestIdx] = min(fit);
        Rabbit_Location = X(bestIdx,:);
        Convergence_curve(t) = Rabbit_Energy;
    end
end

function X = tentInit(N,dim,lb,ub)
    z = rand(N,dim);
    mu = 0.499;
    for k = 1:12
        idx = z < mu;
        z(idx) = z(idx)/mu;
        z(~idx) = (1-z(~idx))/(1-mu);
    end
    X = z.*(ub-lb)+lb;
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

function c = cauchyRand(m,n)
    c = tan(pi*(rand(m,n)-0.5));
    c = max(min(c,20),-20);
end
