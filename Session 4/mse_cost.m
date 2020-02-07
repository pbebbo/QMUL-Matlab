function mse = mse_cost(X, y, w)
mse = mean((X*w-y).^2);
end