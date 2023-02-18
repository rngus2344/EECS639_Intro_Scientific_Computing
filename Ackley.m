function [fx,x1,x2] = Ackley(x1,x2)
%   This function contstructs the Ackley Function for 
%   use in testing unconstrained optimization methods.
%   This function takes inputs of the min and max of the
%   range of x values and outputs strings for x1, x2, and for the 
%   Ackley Function.

    fx = -20*exp(-.2*sqrt(.5*sum([x1^2,x2^2])))-exp(.5*(cos(2*pi*x1)+cos(2*pi*x2)))+20+exp(1);    
end