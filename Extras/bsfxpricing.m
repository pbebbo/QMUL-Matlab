
%Pricing FX options using Black-Scholes
S = 0.8;
K = 0.7;
r = 0.06;
t = 4/12;
Vol = 0.12;
DivRate = 0.08;

[call,put] = blsprice(S,K,r,t,Vol,DivRate);

