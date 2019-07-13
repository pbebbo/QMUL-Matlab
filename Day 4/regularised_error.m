function e = regularised_error(data,coeff,lambda)
e = 0;
for i = 1:length(data)
    px = 0;
    for j = 1:length(coeff)
        px = px + (data(i,1))^(j-1)*coeff(j);
    end
    e = e + 0.5*(px-data(i,2))^2 + 0.5*lambda*coeff'*coeff;
end
end