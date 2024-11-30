%% Matlab Course Session 4

%% Function Handle
fun = @(x)((2*pi).^(-0.5))*exp((-0.5)*x.^2);
figure,
fplot(fun)
set(gca,'fontsize',14)
xlabel('x','fontsize',18)
ylabel('f(x)','fontsize',18,'rotation',360),
title('(({2} {\pi})^{-{0.5}}) {exp}((-{0.5}) {x}^{2})','fontsize',18);


%% Function Strings
figure,
ezplot('(x^2-y^4)*cos(x/y)',[-4*pi,4*pi]),
set(gca,'fontsize',14)
xlabel('x','fontsize',18)
ylabel('y','fontsize',18,'rotation',360)
title('({x}^{2}-{y}^{4}) {cos}({x}/{y}) = {0}','fontsize',18);
hndl=get(gca,'Children');
set(hndl,'LineWidth', 0.5,'Color', [0 0 1]);

%% Portfolio Optimization Example
% Synthetic Data
nmax = 1e4;
[X_synth, xav_synth, NP] = generateSyntheticData(nmax);
[sig_synth, r_synth, w_space_synth] = optimizePortfolio(X_synth, xav_synth, NP);
plotEfficientFrontier(sig_synth, r_synth);

% Real Data
symbols = {'ADBE', 'AAPL', 'MSFT', 'PEP', 'KO', 'AMZN'};
[X_real, xav_real, Nassets] = getRealData(symbols, '01-Jan-2010', datetime('today'));
[sig_real, r_real, w_space_real] = optimizePortfolio(X_real, xav_real, Nassets);
plotEfficientFrontier(sig_real, r_real);

%% Coefficient of Determination
[cor_mat,P_ttest]=corrcoef(X);
coeff_of_determination = cor_mat.^2;

%% Stability of Correlations Matrix
ct1 = zeros(Nassets);
ct2 = zeros(Nassets);
w=0;
for t=1:250:(T-250)
    ct=corrcoef(X(t:t+250,:));
    ct1 = ct1+ct;
    ct2 = ct2+ct.^2;
    w=w+1;
end
P_windows_mean = ct1/w
P_windows_std = sqrt(ct2/w-P_windows_mean.^2)

pp = zeros(Nassets);
for t=1:1000 % permutations
    ct=corrcoef([X(randperm(T),1),X(randperm(T),2),X(randperm(T),3),...
        X(randperm(T),4),X(randperm(T),5),X(randperm(T),6)]);
    pp = pp + (abs(ct)>=abs(cor_mat));
end
P_ttest
P_permutation_test = pp/1000

%% Sharpe Ratio
sharpe_ratio=xav./std(X);
fprintf('sharpe_ratio_asset_1: %.9f\n', sharpe_ratio(1))
fprintf('sharpe_ratio_asset_2: %.9f\n', sharpe_ratio(2))
fprintf('sharpe_ratio_asset_3: %.9f\n', sharpe_ratio(3))
fprintf('sharpe_ratio_asset_4: %.9f\n', sharpe_ratio(4))
fprintf('sharpe_ratio_asset_5: %.9f\n', sharpe_ratio(5))
fprintf('sharpe_ratio_asset_6: %.9f\n', sharpe_ratio(6))


%% MSE Example

X = randn(3,3);
y = randn(3,1);
w = randn(3,1);
mse = mse_cost(X, y, w)
w_star = linreg(X,y)


%% Linear Regression
close all,

X = randn(1000,3);
y = X*[0.3;1.0;-0.3]+randn(1000,1);
mdl = fitlm(X,y)

%% Polyfit Example
clear all, close all,clc
display('Polynomial fitting');
x = [1,2,3,4]';
y = [3,2,0,5]';
figure,
plot(x,y,'*');
%map the data and perform linear regression
for k = 1:4
    xx = basis(x,k);
    w  = linreg(xx,y);
    c  = mse_cost(xx,y,w);
    display(sprintf('Bases dimensions: %g, MSE: %.2f', k, c))
end

%plot the data
xplot = linspace(0,5,100)';
%plot with 2-dimensional basis
k2bases = basis(xplot,2);
xx = basis(x,2);
w2 = linreg(xx,y);
yplot = k2bases*w2;
hold on
plot(xplot,yplot)
%plot with 4-dimenstional basis
k4bases = basis(xplot,4);
xx = basis(x,4);
w4 = linreg(xx,y);
yplot = k4bases*w4;

hold on
plot(xplot,yplot)
hold off

%% Run Unregularised Regression
close, clear, clc;
runRegression(false);

%% Run Regularised Regression
close, clear, clc;
runRegression(true);

%% Sigmoid Function
close, clear, clc;
fplot('1/(1+exp(-u))')


%% Binary Logistic Regression Example
close all; clear all;
warning off;
clc

% Generate synthetic financial data
rng(42); % Set random seed for reproducibility

% Create expanded dataset
N = 100; % Number of samples

% Generate features separately for clarity
yearly_return = randn(N,1)*0.2;          % Yearly return ~ N(0, 0.2)
credit_rating = randi([1,3], N,1);       % Credit rating (1=A, 2=B, 3=C)
debt_equity = rand(N,1)*0.1;             % Debt/Equity
market_cap = rand(N,1)*1000;             % Market cap
current_ratio = rand(N,1)*2 + 1;         % Current ratio

% Combine features into matrix
X = [yearly_return, credit_rating, debt_equity, market_cap, current_ratio];

% Generate target variable (default probability)
prob_default = 1./(1 + exp(-(- yearly_return + ...           % Lower returns increase default
                             (credit_rating-2) + ...         % Lower credit rating increases default
                             debt_equity*10 - ...            % Higher debt/equity increases default
                             (market_cap/1000 - 0.5) - ...   % Lower market cap increases default
                             (current_ratio - 2))));         % Lower current ratio increases default

% Generate binary outcomes
Y = rand(N,1) < prob_default;

% Create dummy variables for credit rating
credit_dummy = dummyvar(credit_rating);
X(:,2) = [];  % Remove original credit rating
X = [X credit_dummy(:,1:end-1)]; % Add dummy variables

% Split data into training and test sets
train_idx = 1:round(0.7*N);
test_idx = (round(0.7*N)+1):N;

X_train = X(train_idx,:);
Y_train = Y(train_idx);
X_test = X(test_idx,:);
Y_test = Y(test_idx);

% Fit logistic regression model
ntrials = ones(length(train_idx),1);
[w, dev, stats] = glmfit(X_train, [Y_train ntrials], 'binomial', 'link', 'logit');

% Make predictions
yfit_train = glmval(w, X_train, 'logit', 'size', ntrials);
yfit_test = glmval(w, X_test, 'logit', 'size', ones(length(test_idx),1));

% Calculate metrics
threshold = 0.5;
y_pred_train = yfit_train > threshold;
y_pred_test = yfit_test > threshold;

% Training metrics
accuracy_train = mean(y_pred_train == Y_train);
precision_train = sum(y_pred_train & Y_train) / sum(y_pred_train);
recall_train = sum(y_pred_train & Y_train) / sum(Y_train);
f1_train = 2 * (precision_train * recall_train) / (precision_train + recall_train);

% Test metrics
accuracy_test = mean(y_pred_test == Y_test);
precision_test = sum(y_pred_test & Y_test) / sum(y_pred_test);
recall_test = sum(y_pred_test & Y_test) / sum(Y_test);
f1_test = 2 * (precision_test * recall_test) / (precision_test + recall_test);

% Display results
fprintf('Training Metrics:\n');
fprintf('Accuracy: %.3f\n', accuracy_train);
fprintf('Precision: %.3f\n', precision_train);
fprintf('Recall: %.3f\n', recall_train);
fprintf('F1 Score: %.3f\n\n', f1_train);

fprintf('Test Metrics:\n');
fprintf('Accuracy: %.3f\n', accuracy_test);
fprintf('Precision: %.3f\n', precision_test);
fprintf('Recall: %.3f\n', recall_test);
fprintf('F1 Score: %.3f\n\n', f1_test);

% Visualizations
figure('Position', [100, 100, 1200, 400]);

% ROC Curve
subplot(1,3,1);
[X_roc,Y_roc,T,AUC] = perfcurve(Y_test, yfit_test, 1);
plot(X_roc, Y_roc);
grid on;
title(sprintf('ROC Curve (AUC = %.3f)', AUC));
xlabel('False Positive Rate');
ylabel('True Positive Rate');

% Feature Importance
subplot(1,3,2);
feature_names = {'Return', 'D/E', 'MCap', 'CRatio', 'CreditA', 'CreditB'};
bar(w(2:end));
xticks(1:length(feature_names));
xticklabels(feature_names);
xtickangle(45);
title('Feature Coefficients');
grid on;

% Prediction Distribution
subplot(1,3,3);
histogram(yfit_test(Y_test==0), 'BinWidth', 0.1, 'Normalization', 'probability');
hold on;
histogram(yfit_test(Y_test==1), 'BinWidth', 0.1, 'Normalization', 'probability');
legend('No Default', 'Default');
title('Prediction Distribution');
xlabel('Predicted Probability');
ylabel('Frequency');
grid on;

% Print model summary
fprintf('Model Summary:\n');
fprintf('Deviance: %.3f\n', dev);
fprintf('Number of observations: %d\n', N);
fprintf('Number of predictors: %d\n', size(X,2));
