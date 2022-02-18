function H0_contour_map(params, levels)
% Compute and plot a contour map of H0 in the theta0-z0 space.
    if nargin < 2
        levels = 40;
    end
    figure
    z0s = linspace(-1,1,1e2+1);
    theta0s = linspace(0,2*pi,1e2);
    [theta0M, z0M] = meshgrid(theta0s, z0s);
    H0s = H_fun(z0M, theta0M, params);
    contour(theta0s, z0s, H0s, levels)
end