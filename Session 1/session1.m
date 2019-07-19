%% Matlab Course Session 1
%  Author Peter A. Bebbington
%
% A course aimmed to teach the basics of using the matlab command window

%% Command Window
why

%% Numerical Commands
-11/3, pi^2, date, 2+3;

%% Variables
(-1)^5;
ans
newvariable = ans^2

%% Arrays
x = [1,2,3,4], y = x*2
x = 1:5; y = x + 1.1
x = 0:10:50; y = x * 0.1

%% Building Matrices
A = [1,2,3; 5,6,7]
B = [1:0.5:3]
M = magic(5)


%% Concatenating
A = [1,2;3,4];
B = [3,4;1,2];
C = [A;B]


%% Matrix Elements Example
a = 1:16;
D = reshape(a,4,4)
D(9)
D(3,1)

%% Special Matrices Example
a = linspace(1,10,5)
length(a), size(a)
b = repmat(a,3,1)

%% Set Function Example
b
unique(b)
c = [1 5.5 10];
a(ismember(a,c))
a
union(a,c)
intersect(a,c)

%% Solving Linear Equations
A = [2 2 -2; 6 4 4; 10 8 6]
v = [10; 2; 8]
e = inv(A)*v

%% Bar Plots
year = [1951:10:2001];
population = [50.23, 52.71, 55.93, 56.93, 57.44, 59.05];
bar(year, population)
axis([1941,2011,40,70])
ylabel('Millions','FontSize',8);
xlabel('Year','FontSize',8);
title('UK population growth','FontSize',14);

%% Line Plots
year = [1951:10:2001];
population = [50.23, 52.71, 55.93, 56.93, 57.44, 59.05];
figure
plot(year,population)
plot(year,population,'kx--')
plot(year,population,'mv-','Markersize',20,'Linewidth',5)
ylabel('Millions','FontSize',14)
xlabel('Year','FontSize',14)
title('UK population growth','FontSize',14)

%% Hold and Subplot
clear
close all
x = [0:0.1:6];
y1 = sin(x); y2 = cos(x);
figure,
plot(x, y1, ':b');
hold on;
plot(x, y2,'--r');hold off;
figure,
subplot(2,1,1); plot(x,y1);
subplot(2,1,2); plot(x,y2);

%% 3D Plots
t = 0:pi/50:10*pi;
figure,
plot3(sin(t),cos(t),t);
xlabel('sin(t)');ylabel('cos(t)');zlabel('t');
grid on
axis square
title('Helix')

%% Mesh, Surface, Contour Plots
x = [0:0.1:4]; y = [0:0.1:4];
[X,Y] = meshgrid(x,y); Z = 2*sin(X) + cos(Y);
figure,
mesh(X,Y,Z),
figure,
surf(X,Y,Z)
figure,
contour(X,Y,Z)

