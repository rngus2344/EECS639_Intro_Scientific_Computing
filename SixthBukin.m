function [fx,x1,x2] = SixthBukin(x1,x2)
%   This function contstructs the Sixth Bukin Function for 
%   use in testing unconstrained optimization methods.
%   This function takes inputs of the min and max of the
%   range of x1 and x2 values and outputs strings for x1, x2, and for the 
%   Sixth Bukin Function.
    
    fx = 100*sqrt(abs(x2-.01*x1^2))+.01*abs(x1+10);  
end