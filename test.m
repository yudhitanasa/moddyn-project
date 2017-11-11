close all, clear all, clc;
syms e f p q u x y z
e= 10^-2; p= 0.5; q=0.05; % condition 1
f=0.1;
u=0 ; % without external excitation
x= -1.20824391947380;
y= 0.580206097986845;
z= -1.20824391947380;
double(subs((x + y - q*x^2 - x*y + u) / e))
double(subs(-y + f*z - x*y + u))
double(subs((x-z) / p))