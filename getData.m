classdef getData<handle
    properties  
        % Parameters about PV
        GC;  %A matrix of power regarding the needs of residents
        GG;  %A matrix regarding the power of PV generation        
        %（Columns correspond to times, rows correspond to users）
    end
    methods  
        function obj=getData(N)  
            dt=2;
 
            load('data.mat'); 
            is_summer=1;
            if is_summer==1
                GC=GC1221*0.001;GG=GG1221*0.001;
            else
                GC=GC621*0.001;GG=GG621*0.001;
            end  
            obj.GC=GC(2:N+1,1:dt:end);% kw
            obj.GG=GG(2:N+1,1:dt:end);% kw  
        end
    end
end