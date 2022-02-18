init;

f = @(h) H0_ode_RHS(h, params);

gammas = linspace(1e-2,1e2,20);
BBars = linspace(1e-2,1-1e-3,20);
uBars = linspace(1e-2,1e2,20);

numhs = 20;

hs = zeros(numhs,length(gammas),length(BBars),length(uBars));
fs = zeros(numhs,length(gammas),length(BBars),length(uBars));

textprogressbar(['Sweeping over parameters: '])
count = 0;
for j = 1 : length(gammas)
    for k = 1 : length(BBars)
        for l = 1 : length(uBars)
            count = count + 1;
            textprogressbar(count / numel(fs) * numhs * 100)
            params.gamma = gammas(j);
            params.BBar = BBars(k);
            params.uBar = uBars(l);
            H0Thresh = g(0,params);
            params.H0Thresh = H0Thresh;
            H0Min = g(pi,params);
            params.H0Min = H0Min;
            hs(:,j,k,l) = linspace(params.H0Min+1e-3, 5*params.H0Thresh,numhs);
            fs(:,j,k,l) = arrayfun(@(h) H0_ode_RHS(h, params), hs(:,j,k,l));
        end
    end
end