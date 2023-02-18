function alpha = BacktrackingLineSearch(fun,x,s)
    %   This function minimizes the step scaling alpha term found in some
    %   optimization methods. It does this by constructing a quadratic function 
    %   q(alpha) in 2-D space, defined by fi(alpha) and alpha, and then
    %   minimizing that function on the specified domain.

    alpha_plot = [0,.5,1];
    
    %   Define 3 Points with the Function and alpha
    for i = 1:length(alpha_plot)
        f(i) = feval(fun,x(1)+alpha_plot(i)*s(1),x(2)+alpha_plot(i)*s(2));
    end
    
    %   Construct Quadratic Function q(alpha)
    p = polyfit(alpha_plot,f,2);
    alpha_plot = linspace(alpha_plot(1),alpha_plot(end),100);
    q = polyval(p,alpha_plot);
    
    %   Find alpha that Minimizes q
    [~,min_ind] = min(q);
    alpha = alpha_plot(min_ind);    
end