%% INTRODUCTION TO MATLAB
%  Session 1 Part 1 --- MATLAB BASICS
% ==============================

% We will work through this file together, following the instructions
% marked "TO DO"

% This is a single line comment, you cannot execute it, but you can use it
% to mark your code (VERY useful when you return to your code after a month

%{
All
these
lines
are also comments, specifically they are
multi-line comments!
%}

%% << Two percentage signs start a code section (see above 'n' below)
disp('Hello World!')
disp('These lines of code are in a code section')
disp('You can run all the code in this section at once.')
disp('Make sure you cursor is in this section')
disp('Then click ''run section'' in the Editor toolbar')

% TO DO 
% Run this code section with <ctrl>+<enter>
% What does <ctrl>+<shift>+<enter> do?

% TO DO 
% run a single line of code by copying it to the clipboard and pasting in the command window
% run a single line of code by selecting it and pressing F9  (<shift>+F7 on Mac)

%% ARRAYS AND DATA TYPES

% clear the memory and all breakpoints (what is a breakpoint?)
clear all;

% The most common data type in Matlab is an array of real numbers.
% To create an empty array, you can use the command zeros.  For example

A = zeros(3,4)

% creates an array full of zeros that has dimensions of 3 x 4.  Note that
% the first index refers to the vertical direction (3), not the horizontal
% direction (4). You can also take a look at it like this: the first
% dimension tells you how many rows there are, and the second how many
% columns.

% There are many other ways to initialize arrays that are full of ones,
% or random numbers

A = 0 * ones(3,4)

%% TO DO - Type the commands into this section, then run the section
% Initialize arrays B, C and D, of sizes (3x3), (4x6) and (2x3x2) using the
% commands 'zeros', 'rand' and 'randn' respectively.  To find out more about
% these commands type 'help zeros' etc. at the prompt

help zeros
help rand
help randn

% create an array E of size (100x200) full of random numbers, using randn

% TO DO
% What happens when you add a semicolon (;) to the end of a line of code?


% TO DO 
% Use the command 'size' to retreive the size of arrays A and B


% TO DO
% Show the contents of A and B by typing their names into the command window


%% Data types (why do we need different data types?)
x = 'q'
xx = uint8( 5 )
z = single( pi )
y = 4.9

% TO DO 
% Use the command 'whos' to investigate which variables are in memory


% Find the 'Workspace' panel in the matlab interface. Verify it contains the same 
% information as the output of 'whos'

% TO DO
% Is the result of the following unexpected? Why is that?
x-y
% Challenge:   How do we fix it? (hint: look up 'double' in the help)

%% TO DO - Run these lines one at a time. 
y1 = 5
% 1y = 5   % What happens here? Why?
y_1 = 5
y = pi
y = i    % Complex number!
y = 1i   % preferred method for complex numbers (easier when viewing)

%% Array creation. Run these lines one at a time
x1 = [1 2 3]
x2 = [1, 2, 3]

x3 = 1:6         % what does the : do?
x4 = 1:2:10      % what do 2 :s do?
x5 = (1:2:10)'   % << what does the apostrope do?

x6 = 10:1:0
x7 = 10:-1:0

% TO DO
% create a vector T containining all the integers between 1 and 1000

% TO DO
% create a vector U containining all the integers between 100 and 500, in reverse order

% TO DO
% create a vector V containining all the multiples of 5 between -100 and 100 inclusive


%% Linspace
% TO DO
% Can you work out what the following commands do?
linspace(2, 10, 30)
linspace(20, 10, 5)
% Look up 'help linspace' if you need!

% TO DO
% Use 'linspace' to generate the following vector 
% [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

% TO DO
% Now use the colon operator : to generate exactly the same sequence


%% Indexing arrays
% Indexing is a way to assign or extract data from i ndividual elements in an 
% array, or a subset of elements in an array.
% Indexing allows us to perform tasks like:
%  - Finding what is the value of the second element of an array
%  - Setting the last element in an array to have the value '20'
%  - Changing all values less than 3 to the value 0
%  - Removing every other element from an array 

%% Indexing 1D arrays
% TO DO 
% - run this section. 
% - Experiment with extracting other values from the array
x = [1:2:10]
x(2)
x(5)

y = [30, 100, 0, -50]
y(4)

%% TO DO - try these lines one at a time. What happens? Why?
x = [100, 200, 300, 400, 500]
% x(0) %
% x(6) %
x(6) = 11 %

%% Indexing ranges of values
x = [100, 200, 300, 400, 500]
x(1:3)
x(3:4)
x(3:end)

x = [100, 200, 300, 400, 500]
x(2:5) = -10  % << Changing multiple values at once

%% TO DO - indexing excersize



%% Strings can also be seen as (i.e. are) 1D arrays
F = 'Hello world'
F(1:5)
F(6:end)

%% Indexing in 2D arrays

% Lets clear all of the data
clear all;
% and create a single array
A = randn(5,7)

% There are 35 elements in this array.  We can index them in two different ways. 

% 1. We can extract elements using full indexing using the format
B = A(2,3)

% TO DO 
% Extract the bottom right element of the matrix A

% TO DO 
% Extract the bottom left element of the matrix A

% TO DO 
% Use 2D indexing to change the value in the second row, fourth column to -30


%% Linear indexing
% Alternatively, we can use linear indexing.  Here we use a single index
% to extract elements.  The elements are counted down the columns first, and
% then across the rows
A
B = A(4)
C = A(15)

% TO DO 
% Set the top right element of the matrix A to be equal to zero using linear
% indexing.


%% More advanced indexing

% Now let's investigate some more complex ways to index arrays
% We can use the colon operator ':' to extract all of the elements in a
% given index position

% TO DO
% Investigate the result of the following commands (uncomment them in turn)
% B = A(:,3)
% C = A(2,:)
% D = A(:)

% We can also extract just a subset of the array

% TO DO 
% Investigate the result of the following commands

% B = A(2:4,:);
% C = A(2:end,:);
% D = A(2:3,4:6);


% We can also use the double colon notation to skip elements

% TO DO
% Investigate the result of the following commands

% B = A(1:3:end,:);
% C = A(2,2:2:6);


% Finally, we can index one array with another
B = zeros(2,1);
B(1) = 4;
B(2) = 8

% TO DO
% Investigate the result of the following commands
% C = A(B)


%% Reshaping arrays
% We can also reshape arrays to different sizes that have the same
% number of elements

A = randn(3,4);

% TO DO 
% Use the command 'reshape' to change the shape of this array to be
% (2x6).  Note the order that the elements wind up in

% TO DO
% Try to reshape A to be (3x3). What happens?


%% Array concatenation (linking in arrays)

% Let's now investigate concatenating arrays
A = randn(2,2)
B = randn(2,2)

% The square brackets notation allows us to concatenate arrays
% A space denotes horizontal concatenation, a semicolon denotes vertical 
% concatentation. 

% TO DO 
% Investigate the result of the following commands

% C = [A B];        % horizontal concatenation
% D = [A;B];        % vertical concatenation
% E = [1 2; 3 4];   % what does this do?
% F = [1 2:2:10; zeros(2,2) randn(2,4)];

% Finally, we can replicate arrays using the command repmat

% TO DO
% Investigate the result of the following commands
% I = repmat(A,3,1)
% J = repmat(B,2,3)


%% Mathematical operators
x = [1:2:10]
y = 23
z = [2:10:42]

% TO DO
% run each of these commands one at a time
% What is happening on each occasion?

% Operators combining a vector (x) and a scalar (y):
A = x+y
B = x-y
C = x/y
D = x*y

% Operators combining a vector (x) with another vector (z)  - i.e. vectorisation
% E = x*z
F = x*z'
G = x'*z
H = x.*z


% TO DO
% What is the difference between the different ways of multiplying two vectors?



%% Congratulations, you can now manage around MATLAB!
% TO DO - move on to exercises


%% INTRODUCTION TO MATLAB
%  Session 1 Part 2 --- MATLAB PLOTTING, FLOW, SETS AND MATRIX OPERATIONS
% ==============================

% Here we will investigate some of the basics of matlab's plotting
% functionality.

%% Line plots
x_values = [0:0.01:2*pi - 0.1]';        % reminder: double colon, apostrophe
sin_values = sin(x_values);             

figure; % Create new one, makes it the CURRENT figure. More about this later.
plot(x_values, sin_values);

%% Changing the plotting style
% This code plots each data point as a blue cross:
plot(x_values, sin_values, 'bx');       % what is 'bx'? how do you find out?

% TO DO
% Plot the same data, but using a red dotted line (hint: use 'help plot')


%% Plotting two datasets on the same graph
cos_values = cos(x_values);

plot(x_values, sin_values, 'b');
hold on                                 % help hold
plot(x_values, cos_values, 'r');
hold off

%% Setting axis limits
% The matlab plot will change the limits of the axes to fit the data
% TO DO Try the following:
tan_values = tan(x_values);

plot(x_values, sin_values, 'b');
hold on
plot(x_values, cos_values, 'r');
plot(x_values, tan_values, 'g');
hold off

% What has happened? Why? (hint - try max(tan_values) to see the maximum value
% of tan_values.

% Set the axis limits manually to something sensible using 'ylim'


%% TO DO - Plot formatting
% Add a title to the plot of all three trigonomic functions.
% Add labels to the x and y axis.
% Add a legend.


%% 3D plotting 1
% Just an example of some more advanced matlab plotting
my3Dfigure = figure; 
surfc(peaks(30))

% TO DO - find the 'rotate 3D' button in the figure window, and view the 3D plot
% from different angles

%% 3D plotting 2
surfl(peaks(30));
lighting phong
shading interp;


%% From here on, we'll talk about operation flow, sets, and matrix operations


%% A basic 'For' loop (repetition for a known number)
for x=1:10
    disp('This statement is repeated')
end

%% Note: this is looping over the COLUMNS of x, so this runs 3 times:
for x = [1 3 5]
    disp('x equals')
    disp(x)
end

% as does this (what is the result?)
for x = [1 2 3; 4 5 6]
    disp('x equals')
    disp(x)
end

%% Performing operation on each element in a loop
A=[1:12]; 
B = reshape(A,4,3)
for x = B
    disp('Mean of column is:'); 
    disp(mean(x)); 
    disp('Sum of column is:'); 
    disp(sum(x));
end

% But note! Many functions already use columns:
sum(B)
mean(B)
std(B)

% TO DO
x = [10, 20, 40, 230];

% Create a for loop to print out the values of x like so:

% Element has value 10
% Element has value 20
% Element has value 40
% Element has value 230

% Hint:
% You can join numbers and strings together like so:
y = 10.3;
message = ['y equals ', num2str(y)];
disp(message)

% TO DO
% Now add some extra elements to x, and check your loop still prints out the
% values

% TO DO
% Now create a for loop to print the following:

% Element number 1 has value 10
% Element number 2 has value 20
% Element number 3 has value 40
% Element number 4 has value 230



%% While loops
% WHILE: based on any expression!
%{
while expression
    statements  % are executed while expression is true
                % if the expression returns an array all
                % must be true
end
%}


a_counter = 0;
while a_counter < 3
    a_counter = a_counter + 1;
    disp('The counter is equal to: ')
    disp(a_counter);
end


%% If statement
help if
%        IF expression
%          statements
%        ELSEIF expression
%          statements
%        ELSE
%          statements
%        END

help switch



%% Basic set operations (set - collection of distinct objects - no duplicates!)

% finding all the unique elements in an array
unique([4 2 2 4 7])

a = [1:5]
b = [7:-1:4]

% set operations and intersection
union(a,b)              % union  (elements which are either in A or B)
intersect(a,b)          % intersection (elements both in A and B)

setdiff(a,b)            % set difference (elements in A, not in B)
setdiff(b,a)            % vice-versa

% TO DO
% Can you find all the elements which are in either a OR b, but not in both,
% using a combination of union, intersect and setdiff?
% Hint: try drawing a Venn Diagram to help you.

%% More complex set operations
ismember(a,b)           % are elements of A in B
a(ismember(a,b))        % reminder: indexing

ismember(b,a)
b(ismember(b,a))


A = [1 2 3 4]; B = [2 3 4 1]; C = [3 4 2 2];

setdiff(A,B),setdiff(B,A)

setdiff(A,C),setdiff(C,A)

A = [1:5],B = [1 2 3 3 4], C = [1:5]

ismember(A,B) == ismember(B,A)

ismember(A,C) == ismember(C,A)

ismember(A,C) == ismember(C,A)



%% Solving matrix operations
% A very powerful capability of matlab is its ability to solve matrix
% operations. Take the following set of equations in three variables:
%
%     2 * x1 + x2 - 2 * x3 = 10
%     6 * x1 + 4 x2 + 4 * x3 = 2
%     10 * x1 + 8 * x2 + 6 * x3 = 8
%
% These can be written as a matrix equation like so:
% Ax = b
%
% Where: (A - coefficients of variables)
A = [2 1 -2;
    6 4 4; 
    10 8 6];

% b - constant terms
b = [10 ; 2 ; 8]

% We can solve for x in the following different ways:
% (remind yourself of matrix inversion!)
%  x = A^-1 * b

% invert A and multiply by b
x = inv(A) * b
x = A^-1 * b

% left division
x = A \ b   
% the difference between left division and the previous methods is 
% that if the matrix is singular, then it will calculate the least squares solution
            
% Using the Moore-Penrose pseudoinverse (also good for singular matrices)
x = pinv(A) * b


% why is this (left division, MP pseudoinverse) important?
% having large, singular matrices, and matrices with small values, you
% will value precision and correctness which inv(A) might not be able to
% achieve

% TO DO
% 
B = [-3 , -5, 1;
    9 ,14, 1 ;
    13, 29, -2]
c = [10; -4; 1]

t1= B\c
t2 = inv(B)*c
t3 = pinv(B)*c

B*t1
B*t2
B*t3

