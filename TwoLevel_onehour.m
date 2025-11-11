classdef TwoLevel_onehour < handle

    properties
        %Decision variable
        PowerDemand;  %[D1,D2,D3]
        PowerSupply;  %[S1,S2,S3]
        
        Powertrana2b; %[]
        Powertrana2c; 
        Powertranb2c;
        PowerGrid;
        
        %Dual variable
        Lambda;       %[lambda1,lambda2,lambda3]
        mu;
        

        %Adjustable parameters
        number_of_Microgrid;
        MA;
        MB;
        MC;

        PVa=1;PVb=1;PVc=1;

        PB=[0,0,0];

        %for mu 
        beta1=0.5;kmaxM=30;
        %for lambda
        beta2=0.5;kmaxL=30;
        %Cost of transmission
        CofT=0.5;

        %show
        %Other intermediate variables
        TotalCost;
    end

    methods
        function obj = TwoLevel_onehour()
            obj.number_of_Microgrid=3;
            obj.MA= Micro_Community_A();
            obj.MB= Micro_Community_B();
            obj.MC= Micro_Community_C();

            obj.Lambda=ones(1,obj.number_of_Microgrid)*0.1; 
        end

        function Solve(obj)
            for k=1:obj.kmaxL 
                obj.UpdatePowerDemand();
 
                obj.UpdatePowerSupply();
 
                lambdaold=obj.Lambda;
                obj.UpdateLambda();

                err1=norm(obj.PowerDemand-obj.PowerSupply+obj.PB);
                err2=norm(obj.Lambda-lambdaold);                
                disp([err1,err2])
                if err1<1e-3 && err2 <1e-3
                    disp(k)
                    break;
                end
            end
            obj.CalculateTotalCost(); 
        end

        function UpdatePowerDemand(obj)
            U =obj.MA.utility_function();
            [DA, ~] = fminbnd(@(x) -U(x)+obj.Lambda(1)*x, 0, obj.MA.PowerDemand_Max);

            U =obj.MB.utility_function();
            [DB, ~] = fminbnd(@(x) -U(x)+obj.Lambda(2)*x, 0, obj.MB.PowerDemand_Max);

            U =obj.MC.utility_function();
            [DC, ~] = fminbnd(@(x) -U(x)+obj.Lambda(3)*x, 0, obj.MC.PowerDemand_Max);

            obj.PowerDemand=[DA,DB,DC];
        end

        function UpdatePowerSupply(obj)
            options = optimoptions('fmincon', 'Display', 'off');
            SA=zeros(1,3);SB=zeros(1,3);SC=zeros(1,3);
            PDMA=obj.MA.PowerDemand_Max;PDMB=obj.MB.PowerDemand_Max;PDMC=obj.MC.PowerDemand_Max;
            mu0=zeros(3,1);

                CgridA =obj.MA.Grid_function(); 
                CgridB =obj.MB.Grid_function(); 
                CgridC =obj.MC.Grid_function();
                Ctran =obj.MA.Tran_function();
            for k=1:obj.kmaxM

                funA = @(x) CgridA(x(1))-Ctran(x(2))-Ctran(x(3))-obj.Lambda(1)*(x(1)+obj.PVa-x(2)-x(3))...
                    +mu0(1)*x(2)+mu0(2)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
                lb=[0,-obj.PVb,-obj.PVc];
                ub=[PDMA,obj.PVa,obj.PVa];
                A=[0,1,1];
                b=obj.PVa;
                [SA, ~] = fmincon(funA, SA, A, b, [], [], lb, ub, [],options);

             
                funB = @(x) CgridB(x(1))+Ctran(x(2))-Ctran(x(3))-obj.Lambda(2)*(x(1)+obj.PVb+x(2)-x(3))...
                    -mu0(1)*x(2)+mu0(3)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
                lb=[0,-obj.PVb,-obj.PVb];
                ub=[PDMB,obj.PVa,obj.PVc];
                A=[0,-1,1];
                b=obj.PVb;
                [SB, ~] = fmincon(funB, SB, A, b, [], [], lb, ub, [],options);

                funC = @(x) CgridC(x(1))+Ctran(x(2))+Ctran(x(3))-obj.Lambda(3)*(x(1)+obj.PVb+x(2)+x(3))...
                    -mu0(2)*x(2)-mu0(3)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
                lb=[0,-obj.PVc,-obj.PVb];
                ub=[PDMC,obj.PVa,obj.PVc];
                A=[0,-1,-1];
                b=obj.PVc;
                [SC, ~] = fmincon(funC, SC, A, b, [], [], lb, ub, [],options);
                
                muold=mu0;
                mu0=mu0+obj.beta1*[SA(2)-SB(2);SA(3)-SC(2);SB(3)-SC(3)];
                
                err1=norm([SA(2)-SB(2);SA(3)-SC(2);SB(3)-SC(3)]);
                err2=norm(muold-mu0);               
                disp(['mu :',num2str([SA(2),SA(3),SB(3)])])
                if err1<1e-3 && err2 <1e-3
                    disp(k)
                    break;
                end
            end

            obj.PowerSupply=[SA(1)+obj.PVa-SA(2)-SA(3),SB(1)+obj.PVb+SB(2)-SB(3),SC(1)+obj.PVc+SC(2)+SC(3)];

            %recode
            obj.PowerGrid=[SA(1),SB(1),SC(1)];
            obj.Powertrana2b=[SA(2),SB(2)];
            obj.Powertrana2c=[SA(3),SC(2)];
            obj.Powertranb2c=[SB(3),SC(3)];
            obj.mu=mu0;
        end
        
        function UpdateLambda(obj)
            obj.Lambda=obj.Lambda+obj.beta2*(obj.PowerDemand-obj.PowerSupply+obj.PB);
        end

        function CalculateTotalCost(obj)
            U1 =obj.MA.utility_function();
            U2 =obj.MB.utility_function();
            U3 =obj.MC.utility_function();

            Cgrid1 =obj.MA.Grid_function();
            Cgrid2 =obj.MB.Grid_function();
            Cgrid3 =obj.MC.Grid_function();

            obj.TotalCost=-(U1(obj.PowerDemand(1))+U2(obj.PowerDemand(2))+U3(obj.PowerDemand(3)));
            obj.TotalCost=obj.TotalCost+(Cgrid1(obj.PowerGrid(1))+Cgrid2(obj.PowerGrid(2))+Cgrid3(obj.PowerGrid(3)));
        end
    end
end