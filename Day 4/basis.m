function phiX = basis(X, k)
phiX = [];
for i=1:k
    phiX = [phiX, X.^(i-1)];
end
end