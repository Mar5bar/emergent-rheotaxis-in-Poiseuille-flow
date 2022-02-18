function phase_portrait(params)
% Compute and plot a phase portrait of the leading-order dynamics.
    figure
    hold on
    ts = linspace(0,10,1e3);
    H0s = linspace(-1,1,1e2);
    theta0Init = pi;
    for H0 = H0s
        try
            z0Init = z_leading_order(theta0Init, H0, params);
            [z0,theta0] = solve_leading_order_odes(z0Init, theta0Init, ts, params);
            plot(mod(theta0,2*pi),z0,'Color','black')
        end
    end
    xlabel('$\theta_0$')
    ylabel('$z_0$')
end