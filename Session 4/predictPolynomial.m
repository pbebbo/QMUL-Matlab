function y_pred = predictPolynomial(x, coeffs)
    % Predict values using polynomial coefficients
    y_pred = zeros(size(x));
    for j = 1:length(coeffs)
        y_pred = y_pred + (x.^(j-1))*coeffs(j);
    end
end 