% QUESTION 7 VERIVY THE POSITIVE INVARIANT FOR UNDRIVEN SYSTEM
close all, clear all, clc;
syms e f p q u x y z
func = @(t,x) [(x(1) + x(2) - q*x(1).^2 - x(1)*x(2) ) / e; -x(2) + f*x(3) - x(1)*x(2) ; (x(1)-x(3)) / p];
[t,xa] = ode45(func,[0 10],[0.5 0.5 0.5]);