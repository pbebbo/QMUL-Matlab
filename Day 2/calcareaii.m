function area = calcareaii(rad)
% This function calculates the area of a circle
% The input argument can be a vector of radii 
% Why are we using .* instead of * ? What happens when rad is a vector?
area = pi * rad .* rad;
end
