%% Fibonacci 
fib = [1 2];% assign intial value to variables
fib_sum = 0;

while ( fib(2) < 4e6 )
    %check if largest value is even
    if ( mod(fib(2),2) == 0)
        fib_sum = fib_sum + fib(2);
    end
    % generate next number in sequence
    fib = [ fib(2) sum(fib) ];
end

%% Exercise 1
% figure out how many terms in the sum 1+2+3+...+n
% exceed one million
array = 1:5000;
size_of_array = numel(array);
sum_of_terms = 0;
number_of_intergers = 0;
while sum_of_terms < 1e6;
    number_of_intergers = number_of_intergers + 1;
    added = array(number_of_intergers);
    sum_of_terms = sum_of_terms + added;
end
number_of_intergers
sum_of_terms

%% Exercise 2
% figure out how many terms in the sum 1+2+3+...+300
% exceed one million
array2 = 1:300;
size_of_array2 = numel(array2);
sum_of_terms = 0;
n = 1;
number_of_intergers2 = 0;
while number_of_intergers2 < size_of_array2;
    number_of_intergers2 = number_of_intergers2 +1;
    added2 = array2(number_of_intergers2);
    sum_of_terms = sum_of_terms + added2;
    if (number_of_intergers2 == n*20);
        n = n + 1;
        number_of_intergers2
        sum_of_terms
    end
end

%% Exercise 3
% write a program to deterime the result of flipping a coin
% 1000 times using a for loop and the rand command
array_coinflip = [];
for coinflip = 1:1000;
    random = rand();
    if random > 0.5;        
        result = 1;
    else random <= 0.5;
        result = 0;
    end
    array_coinflip=[array_coinflip;result];
end

%% Exercise 4
% generate a random walk starting at 1 based on the results of each
% coin flip
array_coinflip = [];
array_coinflip(1) = 1;
result = 1;
for coinflip = 2:1000;
    random = rand();
    if random > 0.5;
        result = result*1.1;
    else random <= 0.5;
        result = result*0.9;
    end
    array_coinflip=[array_coinflip;result];
end
plot(array_coinflip)


