classdef Micro_Community_B

    properties 
        Numbers_of_users=1; 
        PowerDemand_Max=40;
        battery_capacity=20;
    end

    methods
        function obj = Micro_Community_B() 
        end 
        function f=utility_function(obj)
            f = @(x) min(10*x.^0.5,100);
        end
        function f=Grid_function(obj)
            f = @(x) 10*x.^2+1.*x;          
        end
        function f=Tran_function(obj)
            f = @(x) 0.5.*x;          
        end
    end
end