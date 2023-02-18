function type = CritPointClassifier(crit,fun)
    %   This function loops through a list of critical points and
    %   classifies each as a maximum, minimum, saddle, or if there is 
    %   pathological behavior. Pathological points are to be thrown out.
    
    %   Establish Hessian Matrix
    syms x1 x2
    H = hessian(fun,[x1,x2]);
    
    %   Loop Through Crit Points
    for i = 1:length(crit(:,1))
        x = crit(i,3:4);
        x1 = x(1);
        x2 = x(2);
        
        H_k = subs(H);
        
        %   Use Hessian to Determine Type of Critical Point
        if det(H_k) == 0
            type{i} = 'Pathological';
        else
            %   Take the Eigen Values of the Hessian
            check = double(eig(H_k));
            
            %   Use the Eigen Values to Classify the Critical Points
            if all(check > 0) == 1
                type{i} = 'Minimum';
            elseif all(check < 0) == 1
                type{i} = 'Maximum';
            else
                type{i} = 'Saddle';
            end
        end
    end
end