%The number of rolls until you get four '3's.

N = 10000; %Sample size, i.e., number of experiments
M = 200; %Number of throws in each experiment 
p = 1/6; %Probability of getting 3 from rolling a 6-sided die
x = zeros(N,1); 

for i = 1:N
n = 0; %Reset the total throws before getting a 3 for each experiment
x(i) = M; %Initialize the no. of rolls b4 getting a 3 to max no. of throws
    for j = 1:M
        if (rand < p)
            n = n + 1;
            if (n==4) 
                x(i) = j;
            end
        end
    end
end


%Plotting
f = zeros(M,1);
for i = 1:M 
    f(i) = sum(x==i);
end

bar (f/N)
xlabel('Number of rolls');
ylabel('Probability of x rolls');

