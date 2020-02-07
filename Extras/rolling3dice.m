
%Sum of rolling three 6-sided dice

N = 1000; %Sample size
m = 3; %number of dice
s = 6; %Total number of sides of a die
x = zeros(N,1);

dice1 = floor(rand(N,1)*s+1);
dice2 = floor(rand(N,1)*s+1);
dice3 = floor(rand(N,1)*s+1);

x = dice1+dice2+dice3;

%Plotting the PDF
for i = 1:(m*s);
    f(i) = sum(x==i);
end

bar(f/N)
xlabel('Sum of Dice Value');
ylable('Probability');
    

