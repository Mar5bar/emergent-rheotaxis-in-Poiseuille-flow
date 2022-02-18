function bounds = bounds_of_z_oscillations(H0, params)
% Compute and return the bounds of z0 from H0.

    % The max value of z0 is always at theta0 = pi.
    maxBound = z0_fun(pi, H0(:), params);
    % If H0 is below H0Thresh, the minBound is -maxBound.
    mask = H0 < params.H0Thresh;
    minBound = -maxBound .* mask;
    % If H0 is above H0Thresh, the minBound is at theta0 = 0.
    minBound = minBound +  z0_fun(0, H0(:), params) .* ~mask;

    bounds = [maxBound, minBound];
end