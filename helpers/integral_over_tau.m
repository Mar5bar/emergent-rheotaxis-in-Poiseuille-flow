function res = integral_over_tau(H0, params)
% Numerically compute the required integral over a period, as a function of H0.

    % The integration will depend structurally on whether or not H0 is above
    % the threshold H0_thresh. In one case, H0 > H0_thresh, the orbit has z0
    % ~= 0 and goes from theta0 = 0 to theta0 = 2pi. In the other case, the
    % orbit has a restricted range of theta0 and undergoes two time-identical
    % half-loops, with z0 = 0 at each end. Note theta0 is increasing when z0 >
    % 0.
    thresh = params.H0Thresh;

    % We need to compute theta0 and z0 as functions of time.
    if H0 <= thresh
        % We'll go from z0 = 0, theta = thetaMin to z0 = 0, theta = thetaMax.
        thetaMin = gInv(H0, params);
        initCond = [0;thetaMin];
        % Use z0 = 0 as the termination condition.
        eventFun = @lessThanThreshEventFnc;
    else
        % We'll go from theta = 0 to theta = 2*pi.
        thetaMin = 0;
        z0Init = z0_fun(thetaMin, H0, params);
        initCond = [z0Init; thetaMin];
        eventFun = @moreThanThreshEventFnc;
    end

    fun = @(tau, y) [params.uBar * sin(y(2)); params.gamma * y(1) * (1 - params.BBar*cos(2*y(2)))];
    opts = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',eventFun);

    % We don't know how long to integrate for, but we can guess and iterate
    % until we find an event. We could extend the solution, but this won't
    % make too much difference here.
    tauSpan = [0,1/params.gamma];
    sol = struct();
    sol.xe = [];
    while isempty(sol.xe)
        tauSpan(2) = tauSpan(2) * 10;
        sol = ode45(fun, tauSpan, initCond, opts);
    end
    period = sol.xe;
    taus = linspace(0, period, 1e4);
    sol = deval(sol, taus);
    z0 = sol(1, :);
    theta0 = sol(2, :);
    
    % Now compute the required integral.
    integrand = dH0_integrand_dtau(z0, theta0, params);
    res = trapz(taus, integrand);

    % We'll need to double this if H0 <= threshold. Note that the period also
    % doubles.
    if H0 <= thresh
        res = 2 * res;
        period = 2 * period;
    end

    % Finally, divide by the period, as this should be an average.
    res = res / period;

end
% for Hamiltonian theta0. We will compute the appropriate domain for theta0
% and assume that fun(z0,theta0) = - fun(-z0, theta0, params).

% Due to these assumptions, we will take z0 >= 0 always when computing the
% integral, corresponding to the upper branch of the solution (z0, theta0). On
% this branch, theta0 is increasing in time, so we integrate from smallest
% theta0 to biggest theta0. If H0 is larger than the threshold H0Thresh, then
% the trajectory does not loop and we simply integrate. If H0 is smaller than
% the threshold, then we must double our computed integral to account for both
% halves of the loop (noting the antisymmetry of the integrand).
%     thresh = H0_thresh(params);

%     if H0 <= thresh
%         % The range of theta for any H0 less than H0Thresh is given by gInv(H0) and 2*pi - gInv(H0).
%         thetaMin = gInv(H0, params);
%         thetaMax = 2*pi - thetaMin;
%     else
%         % The range of theta is [0,2*pi]
%         thetaMin = 0;
%         thetaMax = 2*pi;
%     end

%     % Integrate over the range of thetas.
%     thetas = linspace(thetaMin,thetaMax,1e4);
%     res = trapz(thetas, fun(z0_fun(thetas, H0, params), thetas, params));

%     if H0 <= thresh
%         % Double the result, as we go over two branches in a loop.
%         % Here, we exploit the oddness of fun in z0.
%         res = 2 * res;
%     end

% end


function [value, isterminal, direction] = lessThanThreshEventFnc(t,y)
    % Terminate when z0 hits 0.
    value = y(1);
    isterminal = 1;
    % z0 will hit 0 from above.
    direction = -1;
end

function [value, isterminal, direction] = moreThanThreshEventFnc(t,y)
    % Terminate when theta0 hits 2*pi.
    value = y(2) - 2*pi;
    isterminal = 1;
    % theta0 will hit 0 from below.
    direction = 1;
end