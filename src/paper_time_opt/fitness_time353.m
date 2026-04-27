function f = fitness_time353(x, prob)
%FITNESS_TIME353 Objective function for time-optimal 3-5-3 trajectory planning.
% x = [t1,t2,t3]. The objective is total running time plus penalties for
% violating joint position, velocity and acceleration constraints.

x = max(min(x(:)', prob.ub), prob.lb);
out = evaluate353(x, prob);

baseCost = sum(x);
penalty = 0;

qViol = max(0, max(out.q - prob.qMax, [], 'all')) + max(0, max(prob.qMin - out.q, [], 'all'));
vViol = max(0, max(abs(out.v) - prob.vMax, [], 'all'));
aViol = max(0, max(abs(out.a) - prob.aMax, [], 'all'));

penalty = penalty + prob.penaltyWeight * (qViol^2 + vViol^2 + aViol^2);

if any(~isfinite(out.q),'all') || any(~isfinite(out.v),'all') || any(~isfinite(out.a),'all')
    penalty = penalty + 1e12;
end

f = baseCost + penalty;
end
