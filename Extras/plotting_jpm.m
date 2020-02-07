%% Plot time series of JPM and corresponding dates in subplot 1
%% Plot ACF/PACF of JPM in subplots 2, 3
figure
subplot(3,1,1)
plot(commonDate,jpm)
datetick('x','mmm-yy')
subplot(3,1,2)
autocorr(jpm,10)
ylim([-0.2 0.2]) % Set limits for y axis
subplot(3,1,3)
parcorr(jpm,10)
ylim([-0.2 0.2]) 

%% ADF test
[h,pVal] = adftest(jpm)

%% Fitting 5 models on JPM and compare
model = {};
model{1} = arima(1,0,0); % AR(1)
model{2} = arima(5,0,0); % AR(5)
model{3} = arima(0,0,1); % MA(1)
model{4} = arima(0,0,5); % MA(5)
model{5} = arima(1,0,1); % ARMA(1,1)

for i = 1:5
    [estModel{i},estParamCov{i},logL{i}] = estimate(model{i},jpm);
end


