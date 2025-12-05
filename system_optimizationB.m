function [G_opt, B_opt,Tab_opt,Tbc_opt] = system_optimizationB(...
    PV, SOC, lambda, sum_x, rho, time_price, ...
    a0, b0, B_max, E_max,tmu,CofT,Tmax)
% 系统优化：求解最优的G和B  
net_demand = sum_x - PV;

options = optimoptions('fmincon', 'Display', 'off', ...
    'Algorithm', 'sqp', 'MaxIterations', 100);


    function f = obje_B(vars)
        G = vars(1);
        B = vars(2);
        Tab=vars(3);
        Tbc=vars(4);

        % 电网成本（考虑分时电价）
        grid_cost = a0*G^2 + (b0 + 0.05*time_price)*G;

        % 电池成本
        if B <= 0  % 充电
            % (1-SOC)/SOC * exp(3) * 1 * (exp(x) - 1)
            battery_cost = (1-SOC)/SOC * exp(0) * 0.1 * (exp(B) - 1);
        else  % 放电
            % 0.01*x^2 + (1-SOC)/SOC * exp(3) * 1 * x
            battery_cost = 0.00011*B^2 + (1-SOC)/SOC * exp(0) * 0.1 * B;
        end 

        % ADMM惩罚项
        admm_penalty = lambda * (-G-PV-B-Tab+Tbc) + ...
            (rho/2) * (sum_x - (G+PV+B+Tab-Tbc))^2;

        %ADMMinner
        admm_inn=tmu(1)*Tab-tmu(3)*Tbc+CofT*Tab^2+CofT*Tbc^2;

        f = grid_cost + battery_cost + admm_penalty+admm_inn;
    end
lb=[0,-min(E_max*(1-SOC),B_max),-Tmax,-Tmax];
ub=[Inf,min(E_max*SOC,B_max),Tmax,Tmax];
x0 = [max(0, net_demand), 0,0,0];
[vars_opt, ~] = fmincon(@obje_B, x0, [], [], [], [], lb, ub, [], options);
G_opt = vars_opt(1);
B_opt = vars_opt(2);
Tab_opt=vars_opt(3);
Tbc_opt=vars_opt(4);  

 
end
 