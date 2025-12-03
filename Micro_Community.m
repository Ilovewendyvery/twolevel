classdef Micro_Community
    properties 
        %Before calculation
        N;T;  
        D_it;
        PV;

        Pmax;
        Gmax=40;
        Tmax=5;
        Bmax=40;
        BC;
        SOC;

        %after calculation
        PowerDemand;PowerSupply;
        Powertran;
        PowerGrid;
        PowerEBSS; 
    end

    methods
        function obj = Micro_Community(data)
            obj.N=size(data.D_it,1);
            obj.T=size(data.D_it,2);
            obj.SOC=zeros(1,obj.T+1);obj.SOC(1)=0.5;

            obj.PV=data.PV;
            obj.D_it=data.D_it;
            obj.Pmax=data.Pmax;

            obj.BC=data.battery_capacity; 
        end 
    end
end
 