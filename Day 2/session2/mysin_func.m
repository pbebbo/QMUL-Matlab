function [y, h] = mysin_func(x)
% mysin takes the input argument and returns the
% sin of the argument. The function also plots
% the x and y values
% calculate the y values
y = sin(x);
% plot the x and y values
h = figure;
plot(x,y);
end