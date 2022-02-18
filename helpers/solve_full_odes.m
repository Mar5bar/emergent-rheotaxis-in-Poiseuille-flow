function [ts, z, theta] = solve_full_odes(zInit, thetaInit, tspan, params)
% Solve the full stiff odes for z and theta.
    f = @(t,z) [params.omega^0.5*params.u(params.omega*t) * sin(z(2)); params.omega^0.5*params.gamma*z(1)*(1 - params.B(params.omega*t)*cos(2*z(2)))];
    opts = odeset('RelTol',1e-9,'AbsTol',1e-9,'OutputFcn',@odetpbar);
    [ts, y] = ode15s(f, tspan, [zInit; thetaInit], opts);
    z = y(:,1);
    theta = y(:,2);
end