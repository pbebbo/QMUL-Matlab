%% Matlab Course Session 2 Part 1


%% Array Operators
X = [10 25 32 47 50 68];
Y = [30 15 8 7 10 28];

addition_X_Y = X + Y;
subtraction_X_Y = X - Y;
sum_multiplication = X * Y';
elementwise_X_Y = X.*Y;
divide_elementwise_X_Y = X./Y;

addition_X_Y
subtraction_X_Y
sum_multiplication
elementwise_X_Y
divide_elementwise_X_Y

%% Boolean Logical Example
a = [3 5 9 0]; b = [2 0 10 1];
(a > 5) & (b <= 10), (a > 5) & (b < 4)
(a > 5) | (b < 4)
xor(a >= 5, b >= 10)
~(a == 3)

%% If-Else Example
a = 57;
if mod(a,7) == 0
    disp('Divisible by 7')
elseif mod(a,7) == 1 || mod(a,7)==6
    disp('Almost divisible by 7')
else
    disp('Definitely not divisible by 7')
end

%% While Example
A = zeros(10,1);
A(7) = 1;
m = 1;
while (A(m) ~= 1)
    m = m + 1;
end

%% Flow Control For

someArray = zeros(100,1);% intialise to optimise code
for i = 1:100
    % the loop runs once for every...
    % elements of the 1:100 array...
    someArray(i) = i+1;
end
% after leaving the the loop i = 100...

%% Data Series
n = 10000;%size of series
series1 = sqrt(2).*randn(n,1);
series2 = sqrt(4).*randn(n,1);

% intialise the variable sum_of_sqr
sum_of_sqr = 0;

% loop through each element
for index = 1:length(series1)
    difference = series1(index) - series2(index);
    sum_of_sqr = sum_of_sqr + difference^2;
end

%% Vectorise
n = 10000;%size of series
series1 = sqrt(2).*randn(n,1);
series2 = sqrt(4).*randn(n,1);

difference = series1 - series2;

sum_of_sqr = difference' * difference;

