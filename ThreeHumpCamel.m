function [fx,x1,x2] = ThreeHumpCamel(x1,x2)
%   This function contstructs the 3 Hump Camel Function for 
%   use in testing unconstrained optimization methods.
%   This function takes inputs of the min and max of the
%   range of x values and outputs strings for x1, x2, and for the 
%   3 Hump Camel Function.
    
    fx = 2*x1^2-1.05*x1^4+x1^6/6+x1*x2+x2^2;
end