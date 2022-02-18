function dH0 = H0_ode_RHS(H0, params)
% Compute and return the RHS of the leading-order ODE for H0.
    
    % We'll need to integrate dH0_integrand over a period in theta0, and
    % multiply the result by -avgIBuMinusUBar.
    dH0 = - integral_over_tau(H0, params) * params.avgIBuMinusUBar;
end