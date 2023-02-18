function [crit,fun_partB] = SteepestDescentMethod()
    %   This function computes the Steepest Descent Method for each of the
    %   given test functions. The backtracking line search method used was
    %   constructing a quadratic function dependent on alpha in 2-D space
    %   and then calculating an alpha to minimize the quadratic function. 
    %   It returns the critical points along with 
    %   their initial guesses on found in the given range of x1 and x2 values.

    %   Specify the Fucntion to be Tested
    fun_test = input('What function to test?\n1) Rosenbrock\n2) 3-hump Camel\n3) Michalewicz\n4) Ackley\n5) 6th Bukin\n6) Test\n');
    clc
    inp = 0;
    tol = 1e-4;
    crit = inf(1,5);
    
    %       Initialize Each Function
    while inp == 0
        switch fun_test
            case 1
                fprintf('Rosenbrock Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = -5:10;
                range_x2 = range_x1;
                fun = @(x1,x2) Rosenbrock(x1,x2);
                fun_partB = 'Rosenbrock';
            case 2
                fprintf('3-hump Camel Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = -5:5;
                range_x2 = range_x1;
                fun = @(x1,x2) ThreeHumpCamel(x1,x2);
                fun_partB = 'Three-Hump Camel';
            case 3
                fprintf('Michalewicz Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = 0:pi;
                range_x2 = range_x1;
                fun = @(x1,x2) Michalewicz(x1,x2);
                fun_partB = 'Michalewicz';
            case 4
                fprintf('Ackley Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = -32.768:32.768;
                range_x2 = range_x1;
                fun = @(x1,x2) Ackley(x1,x2);
                fun_partB = 'Ackley';
            case 5
                fprintf('6th Bukin Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = -15:-5;
                range_x2 = -3:3;
                fun = @(x1,x2) SixthBukin(x1,x2);
                fun_partB = 'Sixth Bukin';
            case 6
                fprintf('Test Function Selected.\nCalculating...\n')
                inp = 1;
                
                range_x1 = -3:3;
                range_x2 = range_x1;
                fun = @(x1,x2) Test(x1,x2);
                fun_partB = 'Test';
                tol = 1e-6;
            otherwise
                fprintf('Please select an available function.\n')
                fun_test = input('What function to test?\n1) Rosenbrock\n2) 3-hump Camel\n3) Michalewicz\n4) Ackley\n5) 6th Bukin\n');
                clc
                inp = 0;
        end
    end
    
    %   Begin Steepest Descent Method       
    %       Iterate to Find Critical Points
    for i = 1:length(range_x1)
        for j = 1:length(range_x2)
            %   Create Initial Guess x0
            x0 = [range_x1(i); range_x2(j)];
            x = x0;
            
%             Calculate Gradient Matrix
            syms x1 x2;
            g = jacobian(fun,[x1, x2]);
            
%             Use x0 to compute Steepest Descent Step
            conv = 1;
            k = 0;
            while conv > tol
%                 Establish x to Evaluate at
                x1 = x(1);
                x2 = x(2);
                
%                 Compute the Step
                g_k = subs(g);
                g_k = g_k';
                s = -double(g_k);
                s = s/norm(s);
                
                %   Utilize Quadratic Function to Minimize alpha to Scale
                %   Step Size
                alpha = BacktrackingLineSearch(fun,x,s);
                
%                 Calculate x at the next step k+1
                x_k1 = double(x+alpha.*s);
                
                if double(isreal(x_k1)) == 0 || norm(double(isnan(x_k1))) ~= 0 || norm(double(isinf(x_k1))) ~= 0
                    break
                end
                
%                 Compute Convergence
                conv = norm(x_k1-x);
                x = x_k1;
                
%                 Update iteration Count
                k = k+1;
            end
            
            %    Record Critical Point
            %        Format = [Initial Guess (x0), Critical Point (x1,x2,f)]
            crit(end+1,:) = [x0(1),x0(2),x(1),x(2),feval(fun,x(1),x(2))];
        end
    end
%     Correct Critical Point Storage Matrix
    crit(1,:) = [];
end
