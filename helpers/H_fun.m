function H0 = H_fun(z, theta, params)
% Return the Hamiltonian corresponding to z and theta with
% params. Note that the leading order approximation the Hamiltonian is just H(z0,theta)
    H0 = params.gamma / (2*params.uBar) * z.^2 + g(theta, params);
end
