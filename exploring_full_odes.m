% Full ODEs for object in Poiseuille flow with varying Bretherton constant.
% dx/dt = V * cos(theta) [decouples]. Note that this is the motion in the frame of the flow, not the lab frame.
% dy/dt = V * sin(theta)
% dtheta/dt = gamma * y * (1 - B * cos(2*theta))

tFinal = 100;
ts = linspace(0,tFinal,1e4);
Ts = linspace(0,1,1e5);

omega = 50;
gamma = 1;
u = @(T) 1 + 0.5*sin(2*pi*T);
B = @(T) 0.32 + 0.3*sin(2*pi*T + 0*pi);
% UB = @(T) U(T).*B(T);

avg = @(fun) trapz(Ts,fun(Ts)) / (max(Ts) - min(Ts));

uBar = avg(u);
BBar = avg(B);


xInit = 0;
yInit = 1;
thetaInit = 0;

f = @(t,z) [omega*u(omega*t) * cos(z(3)); omega*u(omega*t) * sin(z(3)); gamma*z(2)*(1 - B(omega*t)*cos(2*z(3)))];
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
sol = ode15s(f, [0,tFinal], [xInit; yInit; thetaInit], opts);
sol = deval(ts,sol);

x = sol(1,:);
y = sol(2,:);
theta = sol(3,:);

% Solve the leading-order system, just for z0,theta0.
f = @(t,X) [uBar*sin(X(2)); gamma*X(1).*(1-BBar*cos(2*X(2)))];
taus = ts * sqrt(omega);
sol = ode15s(f, [0,taus(end)], [yInit/sqrt(omega); thetaInit], opts);
sol = deval(taus,sol);

z0 = sol(1,:);
theta0 = sol(2,:);

HLeadingOrder = gamma * z0.^2/(2*uBar) + 1 + atanh(sqrt(2*BBar/(1+BBar)) * cos(theta0)) / sqrt(2*BBar*(1+BBar));
HFull = gamma * (y/sqrt(omega)).^2/(2*uBar) + 1 + atanh(sqrt(2*BBar/(1+BBar)) * cos(theta)) / sqrt(2*BBar*(1+BBar));


clf
nexttile()
plot(x,y);
xlabel('$x$')
ylabel('$y$')
% axis equal
nexttile()
plot(ts,y)
xlabel('$t$')
ylabel('$y$')
nexttile()
plot(ts,theta)
xlabel('$t$')
ylabel('$\theta$')
nexttile()
hold on
plot(ts,HLeadingOrder)
plot(ts,HFull)
xlabel('$t$')
ylabel('$H$')
legend({'Leading order','Full numerical sol'})
