function [ts, H0] = solve_H0_ode(H0Init, tspan, params)
% Solve the ODE for H0 from H0Init with params.
    f = @(t,h) H0_ode_RHS(h, params);
    opts = odeset('RelTol',1e-5,'AbsTol',1e-5,'OutputFcn',@odetpbar);
    [ts, H0] = ode45(f, tspan, H0Init, opts);
end