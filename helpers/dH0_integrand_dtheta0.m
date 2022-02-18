function val = dH0_integrand_dtheta0(z0, theta0, params)
% Return the integrand to be used in the ODE for H0 when integrating wrt
% theta0, using the values of z0 and theta0 passed.
    
    val = cos(2*theta0) ./ z0  .* (params.gamma/params.uBar * z0.^2 .* cos(theta0).*(1 - params.BBar*cos(2*theta0)) + sin(theta0).^2);

end