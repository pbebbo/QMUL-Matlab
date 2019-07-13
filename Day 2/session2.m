%% INTRODUCTION TO MATLAB
%  Session 2 Part 1 --- User-defined Functions, Loops and File Handling
% ==============================

% We will work through this file together, following the instructions
% marked "TO DO"

% This is a comment

clc             % reminder: what does this do?
clear all       % and this?

%{
All
these
lines
are comments
%}

%% << Remember, two percentage signs start a code section
disp('Session 2: Introduction to Functions, Loops and File Handling')
disp('We will learn to write user-defined functions, learn how to write efficient code and understand file handling')

% TODO 
% Run this code section with <ctrl>+<enter>
% What does <ctrl>+<shift>+<enter> do?

% TODO 
% run a single line of code by copying it to the clipboard and pasting in
% the command window run a single line of code by selecting it and pressing
% F9  (<shift>+F7 on Mac)

%% PRE-ALLOCATING A VECTOR

% To preallocate a vector is to set the correct size for efficient programming
% Here is an example of how to preallocate a vector

x = zeros(1, 10000)
%x = 0
for i = 2:10000
    x(i) = x(i-1) * 2;
end

% TO DO

% Prints a triangle of stars
% How many will be specified by a variable
% for the number of rows
rows = 3;
for i=1:rows
	for j=1:i           % inner loop just iterates to the value of i (why?)
        fprintf('*')    % what does fprintf do (hint: help)
    end
    fprintf('\n')
end

%% Saving Data to a File
% Prompt the user for rows and columns and
% create a multiplication table to store in
% a file mymulttable.dat

num_rows = input('Enter the number of rows:');
num_cols = input('Enter the number of columns:');
multmatrix = multtable(num_rows, num_cols);
save mymulttable.dat multmatrix -ascii

% TO DO
%  >> createmulttab
% >> load mymulttable.dat
% >> mymulttable


%% TODO - Run these lines one at a time. 
mat = [3:5; 2 5 7]
mymatsum(mat)

% Efficient Method
mat
sum(mat)

% So, to get the overall sum, it is necessary to sum the column sums!

sum(sum(mat))

% TODO take a look at the function outsum = matcolsum(mat)
% matcolsum(mat)

%% Comparisons between efficient and inefficient methods

% % Inefficient method (why?)
% for i = 1:length(v)
%     v(i) = v(i) * 3;
% end
% 
% v
% 
% % Efficient Method (why?)
% v = v*3

%% Logical Operators
% (reminder: ismember)
% For example, let’s say that there is a vector, and we want to compare
% every element in the vector to 5 to determine whether it is greater than
% 5 or not. 
%
%  The result would be a vector (with the same length as the original) with
%  logical true or false values. 

vec = [5 9 3 4 6 11];
isg = vec > 5

%  Notice that this creates a vector consisting of all logical true or
%  false values.
%  Although this is a vector of ones and zeros, and numerical operations
%  can be done on the vector isg, its type is logical rather than double

doubres = isg + 5

% To determine how many of the elements in the vector vec were greater than
% 5, the sum function could be used on the resulting vector isg:
sum(isg) 

% The logical vector isg can also be used to index into the vector. For
% example, if only the elements from the vector that are greater than 5 are
% desired:
vec(isg)

%% Logical Built-in functions
% There are built-in functions in MATLAB that are useful in conjunction
% with vectors or matrices of all logical true or false values; two of
% these are the functions any and all.

% The function any returns logical true if any element in a vector is
% logically true, and false if not.

% The function all returns logical true only if all elements are logically
% true.

% Here are some examples. For the variable vec1, all elements are
% logical true so both any and all return true.

% TO DO
vec1 = [1 3 1 1 2];
any(vec1)               % What is the value? Why?
all(vec1)               % Change vec1(2) = 0. What happens then? Why?

% FIND function 
% The function find returns the indices of a vector that meet some
% criteria. For example, to find all the elements in a vector that are
% greater than 5:
% TO DO
vec = [5 3 6 7 2]
find(vec > 5)

% IS EQUAL Function
% Also, the function isequal is useful in comparing vectors. 
% In MATLAB, using the equality operator with arrays will return 1 or 0 for
% each element; the all function could then be used on the resulting array
% to determine whether all elements were equal or not. The built-in
% function isequal also accomplishes this:

% TO DO
vec1 = [1 3 14 2 99];
vec2 = [1 2 14 3 99];
vec1 == vec2
all(vec1 == vec2)
isequal(vec1,vec2)


%% INTRODUCTION TO MATLAB
%  Session 2 part 2 --- User-defined Functions, Loops and File Handling
% ==============================


%% Vectors and Matrices as Function Arguments
% Using most programming languages, if it is desired to evaluate a function
% on every element in a vector or a matrix, a loop would be necessary to
% accomplish this. However, as we have already seen, in MATLAB an entire
% vector or matrix can be passed as an argument to a function; the function
% will be evaluated on every element. This means that the result will be
% the same size as the argument.

mat = randi([-8 8], 2, 4)           % hint: help
signum(mat)
signmat = signum(mat)
help signum                         % How does this work?

% Vectors or matrices can be passed to user-defined functions as well, as
% long as the operators used in the function are correct. For example, we
% previously defined a function that calculates the area of a circle:

calcareaii(1:3)

% To DO
factgthigh(5000)

% TO DO
% Multiple conditions in a while loop
vec = randn(5,1);
myanywhile(vec)


%% Working with Files

% Using the textscan to read from the file
% TO DO
textscanex

% Using fopen, load and getl functions to read from a file
% TO DO
mat = randi([5 20],2,4)
fid = fopen('randmat.dat','w');
fprintf(fid,'%d %d\n',mat);
fclose(fid);
load randmat.dat

% Writing to Spreadsheet files
% TO DO
ranmat = randi([1 100],5,3)
% xlswrite('ranexcel.xlsx',ranmat)
writematrix(ranmat, 'ranexcel.xlsx')
ssnums = xlsread('ranexcel.xlsx')

%% Using mat files for variables

% To save all variables to a file
save filename
% check using who
who -file filename

% TO DO
mymat =	rand(3,5)
X = 1:6;
Y = X.^2;
who
save examplemat1
who -file examplemat1

% To save just one variable to a file, the format is:
% save filename variablename
%% CELL ARRAYS

% TO DO
% we cannot do this: cellrowvec = [23, 'a', 1:2:9, 'hello']. Why?
cellrowvec = {23, 'a', 1:2:9, 'hello'};

% TO DO
% To create a column vector cell array, the values are instead separated by
% semicolons:
cellcolvec = {23; 'a'; 1:2:9; 'hello'}

% TO DO
% This method creates a 2x2 cell array matrix:
cellmat = {23 'a'; 1:2:9 'hello'}
cellmat{2,1}
cellmat{2,1}(4)

% We can also refer to subsets of cell arrays, for example,
cellcolvec{2:3}
[c1 c2] = cellcolvec{2:3}

% The function cellplot puts a graphical display of the cell array in a
% Figure Window; however, it is a high-level view and basically just
% displays the same information as typing the name of the variable e.g.,
% it wouldn't show the contents of the vector in the previous example

cellplot(cellmat)

%% Structures.. Introduced but fully covered next session

% Create a structure
% TO DO
package = struct('item_no',123,'cost',19.99,...     % what does ... do?
 'price',39.95,'code','g')

% Difference between struct and cell?
package.item_no
package.cost
package.code