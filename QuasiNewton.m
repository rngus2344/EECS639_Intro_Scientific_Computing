function [crit,crit_path] = QuasiNewton(x0_plot)
    %   This function computes the Quasi-Newton Method: BFGS for each of the
    %   given test functions. It returns the critical points along with 
    %   their initial guesses on found in the given range of x1 and x2 values.

    %   Specify the Fucntion to be Tested
    fun_test = input('What function to test?\n1) Rosenbrock\n2) 3-hump Camel\n3) Michalewicz\n4) Ackley\n5) 6th Bukin\n6) Test\n');
    clc
    inp = 0;
    tol = 1e-4;
    crit=inf(1,5);
    
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
                inp =1;
                
                range_x1 = 0:pi;
                range_x2 = range_x1;
                fun = @(x1,x2) Michalewicz(x1,x2);
                fun_partB= 'Michalewicz';
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

    %   Begin BFGS Method       
    %       Iterate to Find Critical Points
    for i = 1:length(range_x1)
        for j = 1:length(range_x2)
            %   Create Initial Guess x0
            x0 = [range_x1(i); range_x2(j)];
            x = x0;
            
            %Create Initial B0
            B0=eye(length(x));
            B = B0;
            
            syms x1 x2  
            %   Calculate Gradient
            g = jacobian(fun,[x1; x2]);
            
            %   Begin Iterations
            conv = 1;
            iter=0;
            while conv > tol
                %   Establish x to Evaluate at
                x1=x(1);
                x2=x(2);
                
                if x1 == -10 && x2 == 1 && fun_partB == "Sixth Bukin"
                       break
                end
                
                % Establishes Gradient of xk
                gOld=double(subs(g))';
                if det(B) < (tol) || B(1,1) < tol
                        break
                end
                
                %Finding sk
                s=-B\gOld;
                
                %Finding xk+1 from sk and xk
                xnew=x+s;
                
                %   Update Convergence Criteria
                conv = norm(xnew-x);
                
                %Establish x to evaluate at
                x1=xnew(1);
                x2=xnew(2);
                
                %Establish gradient of xk+1
                gNew=double(subs(g))';
                
                %Calculate y
                y=gNew-gOld;
                
                %Calculate Bk+1
                B=B+((y*y')/(y'*s))-(B*s*s'*B)/(s'*B*s);
                x=xnew;
                
                %Count iterations
                iter=iter+1;
            end   
            %record critical points
            %        Format = [Initial Guess (x0), Critical Point (x1,x2,f)]
            crit(end+1,:) = [x0(1),x0(2),x(1),x(2),feval(fun,x(1),x(2))];
        end 
    end
    %   Correct Critical Point Storage Matrix
    crit(1,:)=[];

%%   Rerun Code Specifically for Purpose of Plotting Iteration Path
    crit_path = {};
    fun_list = {'Rosenbrock','Three-Hump Camel','Michalewicz','Ackley','Sixth Bukin'};
    
    %       Initialize Each Function    
    %   Run Loop for Each Test Function
    for k = 1:length(fun_list)    
        %   Initialize Function Based on Test Function
        switch fun_list{k}
                case 'Rosenbrock'              
                    range_x1 = -5:10;
                    range_x2 = range_x1;
                    fun = @(x1,x2) Rosenbrock(x1,x2);
                case 'Three-Hump Camel'
                    range_x1 = -5:5;
                    range_x2 = range_x1;
                    fun = @(x1,x2) ThreeHumpCamel(x1,x2);
                case 'Michalewicz'
                    range_x1 = 0:pi;
                    range_x2 = range_x1;
                    fun = @(x1,x2) Michalewicz(x1,x2);
                case 'Ackley'
                    range_x1 = -32.768:32.768;
                    range_x2 = range_x1;
                    fun = @(x1,x2) Ackley(x1,x2);
                case 'Sixth Bukin'
                    range_x1 = -15:-5;
                    range_x2 = -3:3;
                    fun = @(x1,x2) SixthBukin(x1,x2);
                otherwise
                    fprintf('Error. Function name not recognized.\n')
        end
        
%%    Begin BFGS Method             
%%Iterate to Find Critical Points
            %   Create Initial Guess x0
            x0 = x0_plot';
            x = x0;
            
            %Create Initial B0
            B0=eye(length(x));
            B = B0;
            
            syms x1 x2  
            %   Calculate Gradient
            g = jacobian(fun,[x1; x2]);
            
            %   Begin Iterations
            conv = 1;
            iter=0;
            while conv > tol
                %   Establish x to Evaluate at
                x1=x(1,end);
                x2=x(2,end);
                
                if x1 == -10 && x2 == 1 && fun_partB == "Sixth Bukin"
                       break
                end
                
                % Establishes Gradient of xk
                gOld=double(subs(g))';
                if det(B) < (tol) || B(1,1) < tol
                        break
                end
                
                %Finding sk
                s=-B\gOld;
                
                %Finding xk+1 from sk and xk
                xnew=double((x(:,end)+s));
                
                %   Update Convergence Criteria
                conv = norm(xnew-x(:,end));
                
                %Establish x to evaluate at
                x1=xnew(1);
                x2=xnew(2);
                %Establish gradient of xk+1
                gNew=double(subs(g))';
                
                %Calculate y
                y=gNew-gOld;
                
                %Calculate Bk+1
                B=B+((y*y')/(y'*s))-(B*s*s'*B)/(s'*B*s);
                x(:,end+1)=xnew;
                
                %Count iterations
                iter=iter+1;
            end   
            x1_plot{k} = x(1,:)';
            x2_plot{k} = x(2,:)';
            %   Record Critical Points
            crit_path{end+1} = {x1_plot,x2_plot,iter,fun_list{k}};

    end
end
