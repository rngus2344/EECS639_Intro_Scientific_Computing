function [message] = PartBPlotter()
    %   Display Message in Command Window to Show Work in Progress
    fprintf('Plotting in progress...\n')
    
    %   Create 1st Figure of 3-D Plots and Figures of Contour Plots
    %       Create x1 and x2 for Plotting
    fun_list = {'Rosenbrock','Three-Hump Camel','Michalewicz','Ackley','Sixth Bukin'};
    
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
        
        %   Establish Values for Plotting
        x1 = linspace(range_x1(1),range_x1(end),50);
        x2 = linspace(range_x2(1),range_x2(end),50);
        for i = 1:length(x1)
            for j = 1:length(x2)
                f(i,j) = feval(fun,x1(i),x2(j));
            end
        end
    
        %   Create 3-D Plots
        figure(1);
        subplot(3,2,k)
        surf(x1,x2,f)
        xlim([range_x1(1),range_x1(end)])
        ylim([range_x2(1),range_x2(end)])
        title(fun_list{k})
        xlabel('x_1')
        ylabel('x_2')
        zlabel('f(x_1,x_2)')
        set(gca,'FontName','Times New Roman')
        
        figure(k+1)
        contour(x1,x2,f,50);
        xlim([range_x1(1),range_x1(end)])
        ylim([range_x2(1),range_x2(end)])
        xlabel('x_1')
        ylabel('x_2')
        set(gca,'FontName','Times New Roman')
        sgtitle(strcat('Contour Plot of',{' '},fun_list{k},' Function'),'FontName','Times New Roman')
        
        clear f
    end
    figure(1);
    sgtitle('Test Functions','FontName','Times New Roman')
    
    %   Display Message for Command Window to Display Completion
    message = 'Plotting Complete.';
end