%Binomial option pricing
S = 100;
K = 105;
r = 0.045;
t = 15/12;
h = 1/12;
Vol = 0.2;
Type = 1; %Call option
DivRate = 0.02;

[AssetPrice,OptionValue] = binprice(S,K,r,t,h,Vol,Type,DivRate);