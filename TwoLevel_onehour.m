classdef TwoLevel_onehour < handle

    properties
        %Temporary variable
        tSA;tDA; 
        tSB;tDB;
        tSC;tDC;
        tPowerDemand;  %[D1,D2,D3]
        tPowerSupply;  %[S1,S2,S3]
        
        tPowertrana2b; %[a,b]
        tPowertrana2c; %[a,b]
        tPowertranb2c; %[a,b] 
        
        %Dual variable
        lambda;       %[lambda1,lambda2,lambda3] 
        Lambda_convergence=0;
        mu_convergence=0;
        tmu=[0;0;0];  
        mu_ab;%[a,b,c] 
        mu_bc;
        mu_ac;
        

        %Adjustable parameters 
        MA;
        MB;
        MC;
 
 

        %for mu 
        beta1=0.5;kmaxM=50;
        %for lambda
        pp=[0.2,0.2,0.2];kmaxL=500;
        %Cost of transmission
        CofT=0.5;

        %show
        %Other intermediate variables
        TotalCost;
        bool=1;
    end

    methods
        function obj = TwoLevel_onehour(NT) 

            [MA,MB,MC]=date_generation(NT); 
            obj.MA= MA;
            obj.MB= MB;
            obj.MC= MC; 


            obj.lambda=[1,1,1]; 
            obj.tDA=zeros(obj.MA.N,1);
            obj.tDB=zeros(obj.MB.N,1);
            obj.tSA=zeros(1,4);
            obj.tSB=zeros(1,4);
            obj.tSC=zeros(1,4);
        end

        function Solve(obj,t)
            obj.Lambda_convergence=0;
             obj.lambda=[2,2,2]; 
            for k1=1:obj.kmaxL 

                obj.UpdatePowerDemand(t); 

                obj.mu_convergence=0;
                for k2=1:obj.kmaxM
                    obj.UpdatePowerSupply(t);
                    obj.Updatemu(k2);
                    if obj.mu_convergence
                        break;
                    end
                end

                % obj.draw(t);

                obj.UpdateLambda(k1);

                if obj.Lambda_convergence
                    break;
                end 

                
            end
            obj.CalculateTotalCost(t); 
        end

        function UpdatePowerDemand(obj,t) 
            for i=1:obj.MA.N
                b=obj.MA.D_it(i,t);
                [DAi, ~] = fminbnd(@(x) -obj.U_fun(x,b)+obj.lambda(1)*x, 0, b);
                obj.tDA(i)=DAi;
            end
 
            for j=1:obj.MB.N
                b=obj.MB.D_it(j,t);
                [DBj, ~] = fminbnd(@(x) -obj.U_fun(x,b)+obj.lambda(2)*x, 0, b);
                obj.tDB(j)=DBj;
            end 

            b=obj.MC.D_it(1,t);
            [obj.tDC, ~] = fminbnd(@(x) -obj.U_fun(x,b)+obj.lambda(3)*x, 0, b);

            obj.tPowerDemand=[sum(obj.tDA),sum(obj.tDB),sum(obj.tDC)]; 
        end

        function UpdateLambda(obj,k)
            Lambda_old=obj.lambda;
            obj.lambda=obj.lambda+obj.pp.*(obj.tPowerDemand-obj.tPowerSupply); 
            



            err1=norm(obj.tPowerDemand-obj.tPowerSupply);
            err2=norm(obj.lambda-Lambda_old);
            % for kk=1:3
            %     if obj.lambda(kk) < 0.5
            %         obj.pp(kk)=max(0.6*obj.pp(kk),0.00001);
            %     end
            %     if obj.lambda(kk) > 1
            %         obj.pp(kk)=0.2;
            %     end
            % end

            disp([obj.tPowerDemand;obj.tPowerSupply;obj.lambda])
            disp(obj.pp)
            if err1<1e-2 && err2 <1e-2&&k>1
                disp(['iter of lambda====================================:',num2str(k)])
                obj.Lambda_convergence=1;
            end
        end
    
        function UpdatePowerSupply(obj,t)
            options = optimoptions('fmincon', 'Display', 'off');  
            if obj.lambda(1)<0.5
                obj.tSA=[0,0,0,sum(obj.MA.D_it(:,t))-obj.MA.PV(t)];
            else
            
            funA = @(x) obj.G_fun(x(1))-obj.lambda(1)*(x(1)+obj.MA.PV(t)-x(2)-x(3)+x(4))...
                       +obj.E_fun(x(4),obj.MA.SOC(t),-obj.MA.BC*(1-obj.MA.SOC(t)))-obj.tmu(1)*x(2)-obj.tmu(2)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
            lb=[0,-obj.MA.Tmax,-obj.MA.Tmax,-obj.MA.BC*(1-obj.MA.SOC(t))*obj.bool];
            ub=[obj.MA.Gmax,obj.MA.Tmax,obj.MA.Tmax,obj.MA.BC*obj.MA.SOC(t)*obj.bool]; 
            [obj.tSA, ~] = fmincon(funA, obj.tSA, [], [], [], [], lb, ub, [],options);
            end
            
 
            if obj.lambda(2)<0.5
                obj.tSB=[0,0,0,sum(obj.MB.D_it(:,t))-obj.MB.PV(t)];
            else
            funB = @(x) obj.G_fun(x(1))-obj.lambda(2)*(x(1)+obj.MB.PV(t)+x(2)-x(3)+x(4))...
                       +obj.E_fun(x(4),obj.MB.SOC(t),-obj.MB.BC*(1-obj.MB.SOC(t)))+obj.tmu(1)*x(2)-obj.tmu(3)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
            lb=[0,-obj.MB.Tmax,-obj.MB.Tmax,-obj.MB.BC*(1-obj.MB.SOC(t))*obj.bool];
            ub=[obj.MB.Gmax,obj.MB.Tmax,obj.MB.Tmax,obj.MB.BC*obj.MB.SOC(t)*obj.bool];  
            [obj.tSB, ~] = fmincon(funB, obj.tSB, [], [], [], [], lb, ub, [],options);
            end


            if obj.lambda(3)<0.5
                obj.tSC=[0,0,0,sum(obj.MC.D_it(:,t))-obj.MC.PV(t)];
            else
            funC = @(x) obj.G_fun(x(1))-obj.lambda(3)*(x(1)+obj.MC.PV(t)+x(2)+x(3)+x(4))...
                +obj.E_fun(x(4),obj.MC.SOC(t),-obj.MC.BC*(1-obj.MC.SOC(t)))+obj.tmu(2)*x(2)+obj.tmu(3)*x(3)+obj.CofT*(x(2))^2+obj.CofT*(x(3))^2;
            lb=[0,-obj.MC.Tmax,-obj.MC.Tmax,-obj.MC.BC*(1-obj.MC.SOC(t))*obj.bool];
            ub=[obj.MC.Gmax,obj.MC.Tmax,obj.MC.Tmax,obj.MC.BC*obj.MC.SOC(t)*obj.bool];
            [obj.tSC, ~] = fmincon(funC, obj.tSC, [], [], [], [], lb, ub, [],options);
            end


            obj.tPowerSupply=[obj.tSA(1)+obj.MA.PV(t)-obj.tSA(2)-obj.tSA(3)+obj.tSA(4),...
                    obj.tSB(1)+obj.MB.PV(t)+obj.tSB(2)-obj.tSB(3)+obj.tSB(4),...
                    obj.tSC(1)+obj.MC.PV(t)+obj.tSC(2)+obj.tSC(3)+obj.tSC(4)]; 
        end 

        function Updatemu(obj,k)
            mu_old=obj.tmu;
            obj.tmu=obj.tmu-obj.beta1*[obj.tSA(2)-obj.tSB(2);obj.tSA(3)-obj.tSC(2);obj.tSB(3)-obj.tSC(3)]; 

            err1=norm([obj.tSA(2)-obj.tSB(2);obj.tSA(3)-obj.tSC(2);obj.tSB(3)-obj.tSC(3)]);
            err2=norm(obj.tmu-mu_old);
            % disp(['mu err:',num2str([err1,err2])])
            if err1<1e-3 && err2 <1e-3
                disp(['iter of mu:',num2str(k)])
                obj.mu_convergence=1;  
            end
        end



        function CalculateTotalCost(obj,t)
            %recode
            obj.MA.PowerGrid=obj.tSA(1);
            obj.MB.PowerGrid=obj.tSB(1);
            obj.MC.PowerGrid=obj.tSC(1); 
            obj.MA.PowerEBSS=obj.tSA(4);
            obj.MB.PowerEBSS=obj.tSB(4);
            obj.MC.PowerEBSS=obj.tSC(4); 
            obj.mu_ab=obj.tmu(1);obj.mu_ac=obj.tmu(2);obj.mu_bc=obj.tmu(3);

            uA=0;
            for i=1:obj.MA.N
                b=obj.MA.D_it(i,t);
                DAi = obj.U_fun(obj.tDA(i),b);
                uA=uA+DAi;
            end

            uB=0;
            for j=1:obj.MB.N
                b=obj.MB.D_it(j,t);
                DBj=obj.U_fun(obj.tDB(j),b);
                uB=uB+DBj;
            end 

            b=obj.MC.D_it(1,t);
            uC = obj.U_fun(obj.tDC,b);
 

            obj.TotalCost=-(uA+uB+uC);
            obj.TotalCost=obj.TotalCost+(obj.G_fun(obj.tSA(1))+obj.G_fun(obj.tSB(1))+obj.G_fun(obj.tSC(1)));
        end

        function draw(obj,t)
            %画优化图
            clf 
            figure(1)
            subplot(2,2,1)%需求图
            x=linspace(0,2*obj.tDA(1),100);
            hold on;            
            plot(x,obj.U_fun(x,obj.MA.D_it(1,t)),'b-');ymax=max(obj.U_fun(x,obj.MA.D_it(1,t)));
            plot(x,obj.lambda(1)*x,'g-')
            plot([obj.tDA(1), obj.tDA(1)], [0,ymax], 'r-', 'LineWidth', 2); 
            hold off;
            subplot(2,2,2)%供给图
            x=linspace(0,2*obj.tSA(1),100);
            hold on;
            plot(x,obj.G_fun(x),'b-');ymax=max(obj.G_fun(x));
            plot(x,obj.lambda(1)*x,'g-')
            plot([obj.tSA(1), obj.tSA(1)], [0,ymax], 'r-', 'LineWidth', 2);
            hold off;
            subplot(2,2,3)% 
            x=linspace(-2*abs(obj.tSA(4)),2*abs(obj.tSA(4)),100); 
            hold on;
            plot(x,obj.E_fun(x,obj.MA.SOC(t),-obj.MA.BC*(1-obj.MA.SOC(t))));
            ymax=max([obj.E_fun(x,obj.MA.SOC(t),-obj.MA.BC*(1-obj.MA.SOC(t))),obj.lambda(1)*x]);
            ymin=min([obj.E_fun(x,obj.MA.SOC(t),-obj.MA.BC*(1-obj.MA.SOC(t))),obj.lambda(1)*x]);
            plot(x,obj.lambda(1)*x,'g-')
            plot([obj.tSA(4), obj.tSA(4)], [ymin,ymax], 'r-', 'LineWidth', 2);
            hold off;
            disp([obj.tPowerDemand(1),obj.tSA,obj.tPowerDemand(1)-(obj.tSA(1)+obj.MA.PV(t)-obj.tSA(2)-obj.tSA(3)+obj.tSA(4))])
        end
        function draw3(obj,t)
            %画优化图
            clf 
            figure(1)
            subplot(2,2,1)%需求图
            x=linspace(0,2*obj.tDC(1),100);
            hold on;            
            plot(x,obj.U_fun(x,obj.MC.D_it(1,t)),'b-');ymax=max(obj.U_fun(x,obj.MC.D_it(1,t)));
            plot(x,obj.lambda(3)*x,'g-')
            plot([obj.tDC(1), obj.tDC(1)], [0,ymax], 'r-', 'LineWidth', 2); 
            hold off;
            subplot(2,2,2)%供给图
            x=linspace(0,2*obj.tSC(1),100);
            hold on;
            plot(x,obj.G_fun(x),'b-');ymax=max(obj.G_fun(x));
            plot(x,obj.lambda(3)*x,'g-')
            plot([obj.tSC(1), obj.tSC(1)], [0,ymax], 'r-', 'LineWidth', 2);
            hold off;
            subplot(2,2,3)% 
            x=linspace(-2*abs(obj.tSC(4)),2*abs(obj.tSC(4)),100); 
            hold on;
            plot(x,obj.E_fun(x,obj.MC.SOC(t),20));
            ymax=max([obj.E_fun(x,obj.MC.SOC(t),20),obj.lambda(3)*x]);
            ymin=min([obj.E_fun(x,obj.MC.SOC(t),20),obj.lambda(3)*x]);
            plot(x,obj.lambda(3)*x,'g-')
            plot([obj.tSC(4), obj.tSC(4)], [ymin,ymax], 'r-', 'LineWidth', 2);
            hold off;
            disp([obj.tPowerDemand(3),obj.tSC,obj.tPowerDemand(3)-(obj.tSC(1)+obj.MC.PV(t)-obj.tSC(2)-obj.tSC(3)+obj.tSC(4))])
        end
        
    end
    methods (Static)
        function f=U_fun(x,b) 
            beta=100;sigma=1; 
            f = beta*log(1+sigma*x).*(x <= b) + beta*log(1+sigma*b).*(x > b); 
        end
        function f=G_fun(x) 
            f = 1*x.^2+0.02.*x;
        end
        function f=E_fun(x,SOC,nBmax)
            % %SOC=SOC0-x/E;
            alpha=(1-SOC)/SOC*exp(3)*0.1;
            f =  alpha*(exp(x)-1).*(x <= 0) + (0.01*x.^2 + alpha*x).*(x > 0);  
            % f=-1/(2*nBmax)*x.*(x-2*nBmax);
        end  
    end

end
