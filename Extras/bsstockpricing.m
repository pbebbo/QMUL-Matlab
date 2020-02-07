
%Pricing stock options using Black-Scholes
S = 100;
K = 95;
r = 0.1;
t = 0.25;
Vol = 0.5;

[call,put] = blsprice(S,K,r,t,Vol);



