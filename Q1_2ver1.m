%% QUESTION 1 FIND ROOT OF 3 VARIABLES EXPRESSION x, y AND z
close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q=0.05; % condition 1
fm= (0.1:0.1:2); % condition 1
u=0 ; % without external excitation
eq11= (x + y - q*x^2 - x*y + u) / e == 0;% 1st order system x
eq12= -y + f*z - x*y + u == 0; % 1st order system y
eq13= (x-z) / p == 0; % 1st order system z
% system simplification
solzq1= solve(eq13, z);
eq12sub= -y + f*solzq1 - x*y + u == 0; % subtitute z on 1st order system y
solyq1= solve(eq12sub, y);
% then subtitute y on 1st order system x
% the new equation is our system equation
eq11sub= (x + solyq1 - q*x^2 - x*solyq1 + u) / e == 0;
solxq1= solve(eq11sub, x)
% this is the 3 fixed point with 1 fp = o independent of f
% and two fp dependent of f
for i= 1:numel(fm)
    f= fm(i);
    fpq1f(i,1)= double(subs(solxq1(1))); 
    fpq1f(i,2)= double(subs(solxq1(2)));
    fpq1f(i,3)= double(subs(solxq1(3)));
end
ftabq1= horzcat((fm'),fpq1f); % final table with f(x) for question 1
%% QUESTION 2 ASSES STABILITY OF THE FIXED POINTS
% we use eq11sub as our system equation
% to asses stability we have to derive this, then
syms e f p q u x y z % redefine
syms p(x)
p(x)= 100*x - 5*x^2 + (100*f*x)/(x + 1) - (100*f*x^2)/(x + 1);
df= diff(p,x) % this is f'(x)
% for every f value plug all three fp to asses its stability 
fm= (0.1:0.1:2); % condition 1
for i= 1:numel(fm)
    f= fm(i);
    x= fpq1f(i,1); % fixed point 1
    fpq2f(i,1)= double(subs(df));
    x= fpq1f(i,2); % fixed point 2
    fpq2f(i,2)= double(subs(df));
    x= fpq1f(i,3); % fixed point 3
    fpq2f(i,3)= double(subs(df));
end
ftabq2= horzcat(ftabq1,fpq2f); % final table with f'(x) for question 2
