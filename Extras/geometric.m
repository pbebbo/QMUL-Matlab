N = 10000; %Sample size, i.e., number of experiments
M = 100; %Number of throws in each experiment 
p = 1/6; %Probability of getting 6 rolling a 6-sided die
x = zeros(N,1); %Number of times you have to roll a die to get a 6.
for i = 1:N
n = 0; %Reset the total throws before getting a 6 for each experiment
    for j = 1:M
        if (rand < p)
            n = n + 1;
            if (n==1) 
                x(i) = j;
            end
        end
    end
end

sum(x==1)%Number of times you roll the die once to get a 6
sum(x==2)%Number of times you roll the die twice to get a 6

%Plotting, truncating graph at x =25
for i = 1:25 
    f(i) = sum(x==i);
end

bar (f/N)
xlabel('Number of rolls');
ylabel('Probability of x rolls');

