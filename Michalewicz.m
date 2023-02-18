function [fx,x1,x2] = Michalewicz(x1,x2)
%   This function contstructs the Michalewicz Function for 
%   use in testing unconstrained optimization methods.
%   This function takes inputs of the min and max of the
%   range of x values and the k parameter, then it outputs strings for x1, x2, and for the 
%   Michalewicz Function.
    k = 2;
    fx = -(sin(x1)*(sin(1*x1^2/pi))^(2*k))-(sin(x2)*(sin(2*x2^2/pi))^(2*k)); 
end