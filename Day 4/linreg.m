function w = linreg(X,y)
w = (X'*X)\X'*y;
end