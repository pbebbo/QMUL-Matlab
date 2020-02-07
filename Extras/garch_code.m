%% Visual check JPM for cond. heteroskedasticity
figure
subplot(3,1,1)
plot(jpm.^2)
subplot(3,1,2)
autocorr(jpm.^2,10)
subplot(3,1,3)
parcorr(jpm.^2,10)

figure
subplot(3,1,1)
plot(abs(jpm))
subplot(3,1,2)
autocorr(abs(jpm),10)
subplot(3,1,3)
parcorr(abs(jpm),10)


%% Analyze the residuals
% Demean JPM or
% Regress JPM on a constant and obtain the residuals
meanjpm = mean(jpm);
residualjpm = jpm - meanjpm;

% Visual inspection of ARCH effects in the residuals
figure
plot(residualjpm)
title('Residual of JPM')

% Test for ARCH effects in the residuals
% h = 0: no ARCH effect, h = 1: There is ARCH effect
[h,pvalue,stat,cvalue] = archtest(residualjpm,'Lags',1,'Alpha',0.05);

%% Fit various GARCH models onto the real data
% (1) GARCH(p,q)
for i = 1:3
    mdl{i} = garch(i,1);
    [estModel{i},estParamCov{i},logL{i}] = ...
        estimate(mdl{i},residualjpm);
end

% (2) AR(1) + GARCH(1,1) Gaussian innovation
mdl{4} = arima('ARLags',1,'Variance',garch(1,1));
[estModel{4},estParamCov{4},logL{4}] = estimate(mdl{4},jpm);

% (3) AR(1) + GARCH(1,1) Student's t innovation
mdl{5} = arima('ARLags',1,'Variance',garch(1,1));
mdl{5}.Distribution = 't';
[estModel{5},estParamCov{5},logL{5}] = estimate(mdl{5},jpm);


%% Model comparison
numObs = 2721;
for i = 1:5
    numParam(i) = size(estParamCov{i},1);
    [aic(i),bic(i)] = aicbic(logL{i},numParam(i),numObs);
end


%% Infer the conditional variances and residuals
[res,v] = infer(estModel{5},jpm);
% Residuals and cond. variances from fitting Model5 on JPM

figure
subplot(3,1,1)
plot(v)
title('Conditional Variance')

subplot(3,1,2)
plot(res)
title('Residuals')

subplot(3,1,3)
plot(res./sqrt(v))
title('Standardized Residuals')

%% Forecast conditional variances, 10 periods ahead
v_forecast = forecast(estModel{5},10);

%% OTHER STUFFS
%% Simulate time-varying variance and innovation based on GARCH
model = garch('GARCHLags',1,'ARCHLags',1,'Constant',0.1,'GARCH',0.4,'ARCH',0.5);
[v,e] = simulate(model,1000,'NumPaths',1);

% Try to fit a GARCH(1,1) onto the simulated data
modelToEstimate = garch(1,1);
[estModel,estParamCov] = estimate(modelToEstimate,v);

%% Test the significance of autocorrelations
% Ljung-Box Q-test
% H0: All autocorrelations from lags 1-8 are equal to 0
% Try different values for lags
[h,p] = lbqtest(jpm,'Lags',8)
