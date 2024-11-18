%% Matlab Course Session 4

%% Function Handle
fun = @(x)((2*pi).^(-0.5))*exp((-0.5)*x.^2);
figure,
ezplot(fun)
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
set(hndl,'LineWidth',0.5,'Color',[0 0 1]);

%% Portfolio Optimisation - Synthetic Data

nmax = 1e4;% first synthetic data
model1 = garch('Constant',0.01,'GARCH',0.1,'ARCH',0.1);
model2 = garch('Constant',0.01,'GARCH',0.2,'ARCH',0.2);
model3 = garch('Constant',0.01,'GARCH',0.3,'ARCH',0.3);
model4 = garch('Constant',0.01,'GARCH',0.4,'ARCH',0.4);
model5 = garch('Constant',0.01,'GARCH',0.5,'ARCH',0.1);
model6 = garch('Constant',0.01,'GARCH',0.1,'ARCH',0.7);

[x1,rreturns1] = simulate(model1,nmax);
[x2,rreturns2] = simulate(model2,nmax);
[x3,rreturns3] = simulate(model3,nmax);
[x4,rreturns4] = simulate(model4,nmax);
[x5,rreturns5] = simulate(model5,nmax);
[x6,rreturns6] = simulate(model6,nmax);

X  = [rreturns1,rreturns2,rreturns3,rreturns4,rreturns5,rreturns6];
NP = size(X,2);
xav = zeros(1,NP);
for i = 1:NP%normalise signals
    xav(i) = mean(X(:,i));
    X(:,i) = (X(:,i)-xav(i))./std(X(:,i));
end

cov_mat=cov(X);
cor_mat=corrcoef(X);
e = ones(1,NP);

r   = -0.001:0.00001:0.001;
Nr = length(r);
sig = zeros(1,Nr);
b = [zeros(NP,1);1;0];
w_space = zeros(NP,Nr);


A   = [cov_mat e' xav';e 0 0;xav 0 0];

for i=1:Nr                       % generate (r,sig)-space
    b(NP+2)=r(i);
    w=A\b;
    w0=w(1:NP);                   % Optimal wieght
    sig(i)=w0'*cov_mat*w0;        % risk
    w_space(:,i)=w0;
end

figure,
plot(sig,r);
xlabel('Risk, \sigma^2','fontsize',18),
ylabel('Return, R','fontsize',18);

%% Portfolio Optimisation - Real Data

symbols = {'ADBE', 'AAPL', 'MSFT', 'PEP', 'KO', 'NASDX'};
stock_data1 = getMarketDataViaTiingo(symbols, '01-Jan-2010', datetime('today'), 'daily');

Ndays = length(stock_data1(1).AdjClose);
Nassets = length(stock_data1);
X = zeros(Ndays-1,Nassets);
T = Ndays-1;
for i=1:Nassets
    stock_data1(i).Rreturn = price2ret(stock_data1(i).AdjClose);
    X(:,i) = stock_data1(i).Rreturn;
end


xav = zeros(1,Nassets);
for i = 1:Nassets%normalise signals
    xav(i) = mean(X(:,i));
    X(:,i) = (X(:,i)-xav(i))./std(X(:,i));
end

cov_mat=cov(X);
e = ones(1,Nassets);

r   = -0.001:0.0001:0.001;
Nr = length(r);
sig = zeros(1,Nr);
b = [zeros(Nassets,1);1;0];
w_space = zeros(Nassets,Nr);


A   = [cov_mat e' xav';e 0 0;xav 0 0];

for i=1:Nr,                       % generate (r,sig)-space
    b(Nassets+2)=r(i);
    w=A\b;
    w0=w(1:Nassets);                   % Optimal wieght
    sig(i)=w0'*cov_mat*w0;        % risk
    w_space(:,i)=w0;
end

figure,
plot(sig,r);
xlabel('Risk, \sigma^2','fontsize',18),
ylabel('Return, R','fontsize',18);

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
display(sprintf('sharpe_ratio_asset_1: %.9f', sharpe_ratio(1)))
display(sprintf('sharpe_ratio_asset_2: %.9f', sharpe_ratio(2)))
display(sprintf('sharpe_ratio_asset_3: %.9f', sharpe_ratio(3)))
display(sprintf('sharpe_ratio_asset_4: %.9f', sharpe_ratio(4)))
display(sprintf('sharpe_ratio_asset_5: %.9f', sharpe_ratio(5)))
display(sprintf('sharpe_ratio_asset_6: %.9f', sharpe_ratio(6)))


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

%% Unregularised Regression

close all,clear all,clc
%Start of the actual script
noise = 0.05;
%finishing number of data pts
no_data_pts = 100;
%increment in pts
inc = 5;
%Test points
no_test_pts = 100;

%three different polynomials to try and fit with
underfit_coeff = [0,1,1]; %order 2
close_coeff = [0,1,1,1]; %order 3
overfit_coeff = [0,1,1,1,1,1,1,1,1,1]; %order 9

%generate our data to learn
data_pts = generate_data(@f,noise,no_data_pts);
%generate our test data
test_pts = generate_data(@f,noise,no_test_pts);
x = linspace(0,1,1000);
y = f(x);
N_its = length(inc:inc:no_data_pts);
x_RMS = zeros(1,N_its);
y_RMS_data_under = zeros(1,N_its);
y_RMS_data_close = zeros(1,N_its);
y_RMS_data_over  = zeros(1,N_its);
y_RMS_test_under = zeros(1,N_its);
y_RMS_test_close = zeros(1,N_its);
y_RMS_test_over  = zeros(1,N_its);
%Loop through adding more data as we go to get better predictions
for d = inc:inc:no_data_pts
    underfit_coeff = fminsearch(@(c)unregularised_error(...
        data_pts(1:d,:),c),underfit_coeff);
    close_coeff = fminsearch(@(c)unregularised_error(...
        data_pts(1:d,:),c),close_coeff);
    overfit_coeff = fminsearch(@(c)unregularised_error(...
        data_pts(1:d,:),c),overfit_coeff);
    %plot our function predictions
    y_under=zeros(size(x));
    for j = 1:length(underfit_coeff)
        y_under = y_under + (x.^(j-1))*underfit_coeff(j);
    end
    y_close=zeros(size(x));
    for j = 1:length(close_coeff)
        y_close = y_close + (x.^(j-1))*close_coeff(j);
    end
    y_over=zeros(size(x));
    for j = 1:length(overfit_coeff)
        y_over =y_over + (x.^(j-1))*overfit_coeff(j);
    end
    k = d./inc;
    x_RMS(k) = d;
    y_RMS_data_under(k) = RMS_error(data_pts(1:d,:),underfit_coeff);
    y_RMS_data_close(k) = RMS_error(data_pts(1:d,:),close_coeff);
    y_RMS_data_over(k)  = RMS_error(data_pts(1:d,:),overfit_coeff);
    y_RMS_test_under(k) = RMS_error(test_pts,underfit_coeff);
    y_RMS_test_close(k) = RMS_error(test_pts,close_coeff);
    y_RMS_test_over(k)  = RMS_error(test_pts,overfit_coeff);
    
    %figure;
    subplot(2,2,1);
    plot(x,y_under,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=2 Underfitting');
    subplot(2,2,2);
    plot(x,y_close,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=3 close fitting')
    subplot(2,2,3);
    plot(x,y_over,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=9 overfitting')
    subplot(2,2,4);
    plot(x_RMS,y_RMS_data_under,'bx-',x_RMS,y_RMS_test_under,'bx--',...
        x_RMS,y_RMS_data_close,'rx-',x_RMS,y_RMS_test_close,'rx--',...
        x_RMS,y_RMS_data_over,'gx-',x_RMS,y_RMS_test_over,'gx--');
    axis([0,length(data_pts),0,1]);
    title('RMS error')
    legend('d=2 - data', 'd=2 - test','d=3 - data','d=3 - test',...
        'd=9 - data','d=9 - test');
    %print(num2str(d),'-dpng');
    input('Press Enter to continue...\n');
end

%% Regularised Regression
close all,clear all,clc
%Start of the actual script
noise = 0.05;
%finishing number of data pts
no_data_pts = 100;
%increment in pts
inc = 5;
%Test points
no_test_pts = 100;
lambda = exp(-10);
%three different polynomials to try and fit with
underfit_coeff = [0,1,1]'; %order 2
close_coeff = [0,1,1,1]'; %order 3
overfit_coeff = [0,1,1,1,1,1,1,1,1,1]'; %order 9

%generate our data to learn
data_pts = generate_data(@f,noise,no_data_pts);
%generate our test data
test_pts = generate_data(@f,noise,no_test_pts);
x = linspace(0,1,1000);
y = f(x);
N_its = length(inc:inc:no_data_pts);
x_RMS = zeros(1,N_its);
y_RMS_data_under = zeros(1,N_its);
y_RMS_data_close = zeros(1,N_its);
y_RMS_data_over  = zeros(1,N_its);
y_RMS_test_under = zeros(1,N_its);
y_RMS_test_close = zeros(1,N_its);
y_RMS_test_over  = zeros(1,N_its);
%Loop through adding more data as we go to get better predictions
for d = inc:inc:no_data_pts
    underfit_coeff = fminsearch(@(c)regularised_error(...
        data_pts(1:d,:),c,lambda),underfit_coeff);
    close_coeff = fminsearch(@(c)regularised_error(...
        data_pts(1:d,:),c,lambda),close_coeff);
    overfit_coeff = fminsearch(@(c)regularised_error(...
        data_pts(1:d,:),c,lambda),overfit_coeff);
    %plot our function predictions
    y_under=zeros(size(x));
    for j = 1:length(underfit_coeff)
        y_under = y_under + (x.^(j-1))*underfit_coeff(j);
    end
    y_close=zeros(size(x));
    for j = 1:length(close_coeff)
        y_close = y_close + (x.^(j-1))*close_coeff(j);
    end
    y_over=zeros(size(x));
    for j = 1:length(overfit_coeff)
        y_over =y_over + (x.^(j-1))*overfit_coeff(j);
    end
    k = d./inc;
    x_RMS(k) = d;
    y_RMS_data_under(k) = RMS_error(data_pts(1:d,:),underfit_coeff);
    y_RMS_data_close(k) = RMS_error(data_pts(1:d,:),close_coeff);
    y_RMS_data_over(k)  = RMS_error(data_pts(1:d,:),overfit_coeff);
    y_RMS_test_under(k) = RMS_error(test_pts,underfit_coeff);
    y_RMS_test_close(k) = RMS_error(test_pts,close_coeff);
    y_RMS_test_over(k)  = RMS_error(test_pts,overfit_coeff);
    
    %figure;
    subplot(2,2,1);
    plot(x,y_under,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=2 Underfitting');
    subplot(2,2,2);
    plot(x,y_close,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=3 close fitting')
    subplot(2,2,3);
    plot(x,y_over,'b-',x,y,'g-',data_pts(1:d,1),data_pts(1:d,2),'ro');
    axis([0,1,-1.5,1.5]);
    title('d=9 overfitting')
    subplot(2,2,4);
    plot(x_RMS,y_RMS_data_under,'bx-',x_RMS,y_RMS_test_under,'bx--',...
        x_RMS,y_RMS_data_close,'rx-',x_RMS,y_RMS_test_close,'rx--',...
        x_RMS,y_RMS_data_over,'gx-',x_RMS,y_RMS_test_over,'gx--');
    axis([0,length(data_pts),0,1]);
    title('RMS error')
    legend('d=2 - data', 'd=2 - test','d=3 - data','d=3 - test',...
        'd=9 - data','d=9 - test');
    %print(num2str(d),'-dpng');
    input('Press Enter to continue...\n');
end

%% Sigmoid Function
close all,clear all;
clc
ezplot('1/(1+exp(-u))')


%% Binary Logistic Regression:
close all,clear all;
warning off;
clc

X = [0.30 0.05 -0.1 -0.05 0.6  0.0 -0.5  2.0  1.0  -0.4;   %yearly return
     1    2     3    2    1    1    3    2    1     3  ;   %credit rating
     0.02 0.08  0.05 0.09 0.04 0.07 0.02 0.01 0.01  0.03]';%debt/equity   
 
[N,D] =  size(X);
ntrials = ones(N,1);%number on trails
Y = [1 0 0 0 1 0 0 0 0 1]';
credit_dummy= dummyvar(X(:,2));
X(:,2)=[];
X = [X credit_dummy(:,1:end-1)];
[w dev stat] = glmfit( X, [Y ntrials], 'binomial', 'link', 'logit' );
yfit = glmval(w, X, 'logit', 'size', ntrials);
