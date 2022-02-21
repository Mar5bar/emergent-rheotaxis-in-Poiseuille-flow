addpath('./helpers')
% We'll work in the z-theta system.


zInit = 1;
thetaInit = pi/6;

% Range of ts to simulate over.
tspan = [0, 1000];

omega = 50; % Frequency of driving oscillation.
gamma = 1; % Characteristic shear rate.
lambda = 1.2*pi; % Phase shift.
u = @(T) 1 + 0.5*sin(2*pi*T);
B = @(T) 0.32 + 0.3*sin(2*pi*T + lambda);

Ts = linspace(0,1,1e5);
avg = @(fun) trapz(Ts,fun(Ts)) / (max(Ts) - min(Ts));
uBar = avg(u);
BBar = avg(B);

params = struct();
params.zInit = zInit;
params.thetaInit = thetaInit;
params.omega = omega;
params.gamma = gamma;
params.u = u;
params.B = B;
params.uBar = uBar;
params.BBar = BBar;

% We will also need to know the value of an average of IB(T) = int(B(T) - BBar, [0,T]) and (u(T) - uBar).
% We give IB analytically, but one could numerically integrate if needed. Note IB(0) = 0.
IB = @(T) -0.3*cos(2*pi*T + lambda)/(2*pi);
avgIBuMinusUBar = avg(@(T) IB(T).*(u(T) - uBar));
params.avgIBuMinusUBar = avgIBuMinusUBar;
params.avgIuBMinusBBar = -avgIBuMinusUBar;

% The threshold value of H0 for which we swap between tumbling (theta0 in
% [0,2*pi]) and swinging (theta0 in smaller range) is given simply as H0 = g(0,params);
H0Thresh = g(0,params);
params.H0Thresh = H0Thresh;

HInit = H_fun(zInit, thetaInit, params);
params.HInit = HInit;
H0Min = g(pi,params);
params.H0Min = H0Min;

disp(params)
