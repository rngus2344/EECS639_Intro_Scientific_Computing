function crit = NelderMeadMethod(n)
    %   This function computes the Nelder-Mead Optimization Method for a
    %   function of n-dimensions. It returns the critical points along with 
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
    
    %   Establish Initial Parameters
    ReflectParam = 1;
    ExpanParam = 2;
    ContParam = .5;
    ShrinkParam = .5;
    
    for i = 1:length(range_x1)
       for j = 1:length(range_x2)
           %    Create Initial Simplex
           %        Define Initial Vertex of Simplex
           clear x
           x0 = [range_x1(i); range_x2(j)];
           x(:,1) = x0;
           
           if norm(x) >= tol
               h = .05;
           else
               h = 2.5e-4;
           end
           %    Construct the Remaining Vertices of the Simplex
           vertices = 360/(n+1):360/(n+1):360;
           for k = 1:(n+1)
               x(1,k) = x(1,1)+h*cosd(vertices(k));
               x(2,k) = x(2,1)+h*sind(vertices(k));
           end
           
           for k = 1:length(x(1,:))
                   f(k) = feval(fun,x(1,k),x(2,k));
           end
           
           %    Iterate Until Convergence
           conv = 1;
           while conv > tol
               %    1. Re-order x so that f(x1) < f(x2) < ... < f(xn+1)
               count = 0;
               while count < n+1
                   for k = 2:length(f)
                       if f(k) < f(k-1)
                           holder = f(k-1);
                           f(k-1) = f(k);
                           f(k) = holder;

                           holder = x(:,k-1);
                           x(:,k-1) = x(:,k);
                           x(:,k) = holder;                   
                       end
                   end
                   count = count+1;
               end
               
               %    Update Convergence Criteria
               conv = norm([x(:,1);f(1)]-[x(:,end);f(end)]);

               %    2. Calculate Centroid of Simplex (Exclude f(xn+1))
               cx(1,1) = mean(x(1,1:end-1));
               cx(2,1) = mean(x(2,1:end-1));
               fc = feval(fun,cx(1),cx(2));

               %    3. Compute Reflection Step of Transformation
               xr = cx+ReflectParam*(cx-x(:,end));
               fr = feval(fun,xr(1),xr(2));

               %        Insert xr into Simplex if Applicable
               if fr < f(end-1) && fr >= f(1)
                   for k = 2:length(f(1:end-1))
                       if fr < f(k)
                           f((k+1):end) = f(k:(end-1));
                           x(:,(k+1):end) = x(:,k:(end-1));

                           f(k) = fr;
                           x(:,k) = xr;
                       end
                   end
               %    4. Expansion Step of Transformation    
               elseif fr < f(1)
                   xe = cx+ExpanParam*(xr-cx);
                   fe = feval(fun,xe(1),xe(2));

                   %        Insert xe into Simplex if Applicable
                   f(end) = fe;
                   x(:,end) = xe;

               %    5. Contraction Step of Transformation
               elseif fr >= f(end-1)
                   xcon = cx+ContParam*(x(:,end)-cx);
                   fcon = feval(fun,xcon(1),xcon(2));

                   %    Insert xcon into Simplex if Applicable
                   if fcon < f(end)
                       f(end) = fcon;
                       x(:,end) = xcon;
                       
                   %    6. Shrink Contraction Step of Transformation
                   else
                       for k = 2:length(x(1,1:end))
                           x(:,k) = x(:,1)+ShrinkParam*(x(:,k)-x(:,1));
                           f(k) = feval(fun,x(1,k),x(2,k));
                       end                     
                   end
               end
           end
           
           %    Record Critical Point
           %        Format = [Initial Guess (x0), Critical Point (x1,x2,f)]
           crit(end+1,:) = [x0(1),x0(2),x(1,1),x(2,1),f(1)];
       end
    end
    crit(1,:) = [];
end