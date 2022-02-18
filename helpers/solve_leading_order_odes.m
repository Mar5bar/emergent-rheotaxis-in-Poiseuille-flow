function [z0, theta0] = solve_leading_order_odes(z0Init, theta0Init, ts, params)
% Solve the leading-order system of ODEs numerically from given initial
% conditions, sampling at the times in ts, where ts(1) = 0.
    f = @(t,y) [params.uBar * sin(y(2)); params.gamma * y(1) * (1 - params.BBar*cos(2*y(2)))];
    opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
    sol = ode45(f, ts, [z0Init; theta0Init], opts);
    sol = deval(ts,sol);
    z0 = sol(1,:);
    theta0 = sol(2,:);
end