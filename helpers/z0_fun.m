function z0 = z0_fun(theta0, H0, params)
% Return z0 as computed from theta0 and the leading-order Hamiltonian H0. We
% will take the positive square root here, but note that the negative root may
% be required in practice. Numerical errors may result in imag(z0)~=0, so we
% take the real part to eliminate this easily.
    z0 = real(sqrt(2*params.uBar / params.gamma * (H0 - g(theta0, params))));
end