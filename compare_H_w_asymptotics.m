disp("Make sure problem has been set up in init.m")
% Run init.
init

% We now have all the information we need to run a comparison.

% Compute the asymptotic approximation to H0 via the scalar ODE.
disp("Solving asymptotic problem:")
[tsAsym, HAsym] = solve_H0_ode(HInit, tspan, params);

% Compute the solution to the full z-theta system.
disp("Solving full problem:")
[tsFull, zFull, thetaFull] = solve_full_odes(zInit, thetaInit, tspan, params);
% Compute H from this full solution.
HFull = H_fun(zFull, thetaFull, params);

% Plot the full Hamiltonian and the approximation on the same axes;
figure
hold on
plot(tsFull, HFull)
plot(tsAsym, HAsym,'LineWidth',2,'Color','black')

xlabel('$t$')
ylabel('$H$')
legend({'Full solution','Leading-order approx'})

% Plot the evolution of z over time, predicting the amplitude of eventual
% oscillations from the asymptotic result.
zBoundsAsym = bounds_of_z_oscillations(HAsym, params);

figure
hold on
handles = [];
handles = [handles; plot(tsFull, zFull)];
handles = [handles; plot(tsAsym, zBoundsAsym,'LineWidth',2,'Color','black')];
xlabel('$t$')
ylabel('$z$')
legend(handles(1:2),{'Full solution','Asymptotic bounds prediction'},'location','best')