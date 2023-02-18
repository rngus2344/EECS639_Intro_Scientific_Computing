function [] = FindGlobalMinimum()
    %   Specify the Fucntion to be Tested
    fun_test = input('What function to test?\n1) Michalewicz\n2) Ackley\n3) 6th Bukin\n');
    clc
    inp = 0;
    
    %       Initialize Each Function
    while inp == 0
        switch fun_test
            case 1
                fprintf('Michalewicz Selected.\nCalculating...\n')
                inp = 1;
                fun_partB = 'Michalewicz';
                
                x0 = [1,1];
                fun = @(x0) -sin(x0(1)^2/pi)^4*sin(x0(1)) - sin((2*x0(2)^2)/pi)^4*sin(x0(2));
            case 2
                fprintf('Ackley Selected.\nCalculating...\n')
                inp = 1;
                fun_partB = 'Ackley';
                
                x0 = [1,1];
                fun = @(x0) -20*exp(-.2*sqrt(.5*sum([x0(1)^2,x0(2)^2])))-exp(.5*(cos(2*pi*x0(1))+cos(2*pi*x0(2))))+20+exp(1);
            case 3
                fprintf('6th Bukin Selected.\nCalculating...\n')
                inp = 1;
                fun_partB = 'Sixth Bukin';
                
                x0 = [1,1];
                fun = @(x0) 100*sqrt(abs(x0(2)-.01*x0(1)^2))+.01*abs(x0(1)+10);
            otherwise
                fprintf('Please select an available function.\n')
                fun_test = input('What function to test?\n1) Rosenbrock\n2) 3-hump Camel\n3) Michalewicz\n4) Ackley\n5) 6th Bukin\n');
                clc
                inp = 0;
        end
    end
    
%     Find global minimum
    [x,fval] = fminunc(fun,x0)
end