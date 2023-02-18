function [fx,x1,x2] = Rosenbrock(x1,x2)
%   This function contstructs the Rosenbrock Function for 
%   use in testing unconstrained optimization methods.
%   This function takes inputs of the min and max of the
%   range of x values and outputs strings for x1, x2, and for the 
%   Rosenbrock Function.
    

fx = 100*(x2-x1^2)^2+(x1-1)^2;
    
end