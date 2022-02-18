function res = gInv(val, params)
% Compute the inverse of g applied to val, giving a value of theta in [0,pi].
    res = acos(tanh(sqrt(2*params.BBar*(1+params.BBar)) * val) / sqrt(2*params.BBar / (1+params.BBar)));
end