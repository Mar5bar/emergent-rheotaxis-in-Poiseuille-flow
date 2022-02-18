function phase_portrait(params)
% Compute and plot a phase portrait of the leading-order dynamics.
    figure
    hold on
    ts = linspace(0,100,1e3);
    H0s = linspace(params.H0Min,2*params.H0Thresh,1e1);
    theta0Init = pi;
    for H0 = H0s
        for sig = [-1,1]
            z0Init = sig * z0_fun(theta0Init, H0, params);
            [z0,theta0] = solve_leading_order_odes(z0Init, theta0Init, ts, params);
            i = 1;
            while i <= length(theta0)
                if theta0(i) >= 2*pi
                    theta0(i:end) = theta0(i:end) - 2*pi;
                    theta0 = [theta0(1:i-1), NaN, theta0(i:end)];
                    z0 = [z0(1:i-1), NaN, z0(i:end)];
                end
                if theta0(i) < 0
                    theta0(i:end) = theta0(i:end) + 2*pi;
                    theta0 = [theta0(1:i-1), NaN, theta0(i:end)];
                    z0 = [z0(1:i-1), NaN, z0(i:end)];
                end
                i = i + 1;
            end
            plot(theta0,z0,'Color','black')
        end
    end
    xlabel('$\theta_0$')
    ylabel('$z_0$')
    xlim([0,2*pi])
end