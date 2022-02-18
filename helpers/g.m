function res = g(theta, params)
    res = atanh(sqrt(2*params.BBar/(1+params.BBar)) * cos(theta)) / sqrt(2*params.BBar*(1+params.BBar));
end