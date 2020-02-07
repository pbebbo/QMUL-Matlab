
%%% Pricing American put option using Binomial method %%%
%%% Option specifications %%%
S = 45;
K = 45;
r = 0.05;
sigma = 0.2;
T = 3;

%%% Binomial parameters %%%
N = 512; %Number of interval
dt = T/N; %Size of time step
A = 0.5*(exp(-r*dt) + exp((r + sigma^2)*dt));
up = A + sqrt(A^2 -1);
down = 1/up;
p = (exp(r*dt) -down)/(up-down);

%%% Branching the tree %%%
BinTree = S*up.^((N:-1:0)').*down.^((0:N)');

%%% Compute option value at maturity %%%
BinTree = max(K-BinTree,0);

%%% Compute expectations %%%
Expect = p*eye(N+1,N+1) + (1-p)*diag(ones(N,1),1);

for i = N:-1:1 %Going back one step at a time
    DiscountedVal = Expect(1:i,1:i+1)*BinTree*exp(-r*dt);
    %Discounting value one step back
    shareprices = S*up.^((i-1:-1:0)').*down.^((0:i-1)');
    %Share prices at time step
    BinTree = max(DiscountedVal, K-shareprices);
    %Comparing discounted vs intrinsic value for exercise decision
end

Option_price = BinTree;
display('The option price is'),display(Option_price)

