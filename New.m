classdef New < handle
    properties
        %% 参数设置
        n;           % 每个微网用户数量
        T = 24;          % 时间周期数（24小时）
        beta = 15;      % 效用函数系数
        sigma = 8;     % 效用函数形状参数

        % 电网成本参数
        a0 = 0.01;       % 二次项系数 (¥/kW²)
        b0 = 0.2;        % 一次项系数 (¥/kW)
        c0 = 0;          % 常数项

        % 电池参数
        E_max = 1000;     % 电池容量 (kWh)
        B_max = 20;      % 最大充放电功率 (kW)
        SOC_min = 0.15;  % 最小SOC
        SOC_max = 0.95;  % 最大SOC
        SOC_initial = 0.5; % 初始SOC

        %电车参数
        E_ev=100;

        %传递参数
        CofT=0.05;
        Tmax=10;

        % ADMM参数
        rho = 0.8;       % 初始惩罚参数
        max_admm_iter = 200; % 每个时段最大ADMM迭代次数
        abs_tol = 1e-5;
        rel_tol = 1e-5;


        PV;        %(3,T)
        D1;D2;D3;  % (n,T)      

        X1_opt;      % (n,T)最优用电量 [kW]
        X2_opt;      % (n,T)最优用电量 [kW]
        X3_opt;      % (n,T)最优用电量 [kW]
        G_opt;      % (3,T)最优电网供电 [kW]
        B_opt;      % (3,T)最优电池充放电 [kW]（正为放电，负为充电）
        SOC_opt;    % (3,T)SOC轨迹
        Lambda_opt; % (3,T)最优对偶变量
        Tab_opt;
        Tac_opt;
        Tbc_opt;

        total_cost_all;
        covnerge_numer;%(T,1)
        converge_error_lam;%(T,3)
        converge_error_mu;%(T,3)
    end

    methods
        function obj = New()
            %生成模拟数据
            obj.T=24; obj.n=[60,1,1];
            rng(2023); % 设置随机种子 

            %% For A B C
            [D,PV]=generationA(obj.n(1),obj.T);obj.D1=D;obj.PV(1,:)=PV;
            [D,PV]=generationB(obj.n(2),obj.T);obj.D2=D;obj.PV(2,:)=PV;
            [D,PV]=generationC(obj.n(3),obj.T);obj.D3=D;obj.PV(3,:)=PV;

        end
        function run_inital(obj)
            %% 存储结果
            obj.X1_opt = zeros(obj.n(1), obj.T);      % 最优用电量 [kW]
            obj.X2_opt = zeros(obj.n(2), obj.T);      % 最优用电量 [kW]
            obj.X3_opt = zeros(obj.n(3), obj.T);      % 最优用电量 [kW]
            obj.G_opt = zeros(3, obj.T);      % 最优电网供电 [kW]
            obj.B_opt = zeros(3, obj.T);      % 最优电池充放电 [kW]（正为放电，负为充电）
            obj.SOC_opt = zeros(3, obj.T);    % SOC轨迹
            obj.Lambda_opt = zeros(3, obj.T); % 最优对偶变量

            obj.Tab_opt = zeros(1, obj.T);
            obj.Tac_opt = zeros(1, obj.T);
            obj.Tbc_opt = zeros(1, obj.T);

            obj.covnerge_numer = zeros(obj.T, 1); % 存储每个时段的收敛信息  
            obj.total_cost_all = 0;
        end

        function run(obj)
            tic;
            run_inital(obj); 
            SOC_current = obj.SOC_initial*ones(3,1);
            SOC_ev = obj.SOC_initial*ones(obj.n(1)/2,1);
            %% 主循环：按时间顺序处理（时间层在外层）
            for t = 1:obj.T
                obj.hour(t,SOC_current,SOC_ev);

                % 更新SOC（供下一个时段使用）
                [SOC_current,SOC_ev] = update_battery_SOC(SOC_current, obj.B_opt(:,t), ...
                    obj.E_max, obj.SOC_min, obj.SOC_max,SOC_ev,obj.X1_opt(obj.n(1)/2+1:end, t),obj.E_ev);
            end 
            toc; 
        end
        
        function hour(obj,t,SOC_current,SOC_ev)
            fprintf('\n=== 时段 %02d:00 - %02d:00 ===\n', t-1, t);
            fprintf('光伏预测: %.2f kW, 分时电价: %.2f ¥/kWh\n', obj.PV(1,t));
            

            %% ADMM初始化（当前时段）
            rho_t = obj.rho*ones(3,1);  % 当前时段的惩罚参数
            x1 = obj.D1(:, t);
            x2 = obj.D2(:, t);
            x3 = obj.D3(:, t);      % 初始化为基准需求
            G = [max(0, sum(x1) - obj.PV(1,t));max(0, sum(x2) - obj.PV(2,t));max(0, sum(x3) - obj.PV(3,t));]; % 初始电网供电
            B = [0;0;0];        % 初始电池功率
            lambda = [0;0;0];   % 对偶变量
            z_prev=G+obj.PV(:,t)+B;

            %% 存储迭代历史（用于分析和调试）
            primal_res_history = zeros(3,obj.max_admm_iter);
            dual_res_history = zeros(3,obj.max_admm_iter);
            cost_history = zeros(3,obj.max_admm_iter);

            %% ADMM迭代（当前时段）
            converged = false;
            iter_used = obj.max_admm_iter;

            for k = 1:obj.max_admm_iter

                % 保存上一次的z值（用于计算对偶残差）
                if k > 1
                    z_prev = z;
                end

                %% 步骤1: 更新x（用户侧，可并行）
                x1=argmin_x(obj,x1,G(1),B(1),lambda(1), rho_t(1),t,[obj.n(1)/2,obj.n(1)/2],obj.D1(:, t),SOC_ev);
                x2=argmin_x(obj,x2,G(2),B(2),lambda(2), rho_t(2),t,[obj.n(2),0],obj.D2(:, t),[]);
                x3=argmin_x(obj,x3,G(3),B(3),lambda(3), rho_t(3),t,[obj.n(3),0],obj.D3(:, t),[]);

                 

                %% 步骤2: 更新系统变量 (z, G, B)
                [G,B,z,Tab,Tac,Tbc,conver_mu]=argmin_other2(obj,x1,x2,x3,lambda, rho_t,t,SOC_current);



                %% 步骤3: 更新对偶变量
                primal_res = [sum(x1);sum(x2);sum(x3)] - z;
                lambda = lambda + rho_t .* primal_res;

                %% 计算收敛指标
                primal_res_history(:,k) = abs(primal_res);

                if k > 1
                    dual_res_history(:,k) = rho_t .* abs(z - z_prev);
                else
                    dual_res_history(:,k) = Inf*[1;1;1];
                end

                % 计算当前解的总成本（用于监控）
                [u_utility,g_cost,b_cost]=obj.cost_compute(x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current);
                cost_history(:,k) = -u_utility + g_cost + b_cost; 

                %% 收敛检查
                eps_primal = obj.abs_tol + obj.rel_tol * max(norm([sum(x1);sum(x2);sum(x3)]), norm(z));
                eps_dual = obj.abs_tol + obj.rel_tol * norm(z);

                if k > 1 && norm(primal_res_history(:,k)) < eps_primal && ...
                        norm(dual_res_history(:,k)) < eps_dual
                    fprintf('  ADMM收敛于 %d 次迭代\n', k);
                    fprintf('  原始残差: %.2e, 对偶残差: %.2e\n', ...
                        norm(primal_res_history(:,k)), norm(dual_res_history(:,k)));
                    converged = true;
                    iter_used = k;
                    break;
                end

                %% 自适应调整rho（每10次迭代）
                if mod(k, 10) == 0 && k > 10
                    mu = 10;  % 调整阈值
                    tau = 1.5; % 调整因子
                    for kk=1:3
                    if primal_res_history(kk,k) > mu * dual_res_history(kk,k)
                        rho_t(kk) = rho_t(kk) * tau;
                        fprintf('  迭代 %d: 增加rho到 %.3f\n', k, rho_t);
                    elseif dual_res_history(kk,k) > mu * primal_res_history(kk,k)
                        rho_t(kk) = rho_t(kk) / tau;
                        fprintf('  迭代 %d: 减少rho到 %.3f\n', k, rho_t);
                    end
                    end
                end

                %% 显示进度
                if mod(k, 20) == 0 && k < obj.max_admm_iter
                    fprintf('  迭代 %3d: 原始残差=%.2e, 对偶残差=%.2e, 成本=%.2f\n', ...
                        k, norm(primal_res_history(:,k)), norm(dual_res_history(:,k)), norm(cost_history(k)));
                end
            end

            if ~converged
                fprintf('  ADMM未收敛，使用最后迭代结果\n');
            end

            %% 存储当前时段最优解
            obj.X1_opt(:, t) = x1;obj.X2_opt(:, t) = x2;obj.X3_opt(:, t) = x3;
            obj.G_opt(:,t) = G;
            obj.B_opt(:,t) = B;
            obj.SOC_opt(:,t) = SOC_current;
            obj.Lambda_opt(:,t) = lambda;
            obj.Tab_opt(t)=Tab;
            obj.Tac_opt(t)=Tac;
            obj.Tbc_opt(t)=Tbc;

            % 记录收敛信息
            obj.covnerge_numer(t) = iter_used;
            if t==3
            obj.converge_error_lam=[sqrt(sum(primal_res_history(:,2:iter_used).^2,1)); sqrt(sum(dual_res_history(:,2:iter_used).^2,1))];
            obj.converge_error_mu=conver_mu;
            end
            % obj.covnerge_numer(t, 2) = max(primal_res_history(:,min(iter_used, obj.max_admm_iter)));

            % 计算当前时段成本

            [u_utility,g_cost,b_cost]=obj.cost_compute(x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current);
            cost_history(:,k) = -u_utility + g_cost + b_cost;

            obj.total_cost_all = obj.total_cost_all + cost_history(:,k);

            fprintf('  优化结果: 总需求=%.2f kW, 电网=%.2f kW, 电池=%.2f kW, SOC=%.3f\n', ...
                sum(x1), G(1), B(1), SOC_current(1));
            fprintf('  时段成本: 用户效用=%.2f, 电网成本=%.2f, 电池成本=%.2f, 总成本=%.2f\n', ...
                u_utility(1), g_cost(1), b_cost(1), obj.total_cost_all(1)); 

        end

        function x=argmin_x(obj,x,G,B,lambda, rho_t,t,n,D,SOC_ev)            
            % 并行更新每个用户的用电量 
            for i = 1:n(1)
                % 用户i的本地优化问题
                z0=sum(x)-x(i)-(G+obj.PV(t)+B);
                fun = @(xi) user_local_optimization(xi, D(i), lambda, z0, ...
                    obj.beta, obj.sigma, rho_t,'other');

                % 边界约束
                lb = 0;
                ub = 1 * D(i); % 最大不超过2.5倍基准需求

                % 求解一维优化问题
                options = optimset('Display', 'off', 'TolX', 1e-6);
                x(i) = fminbnd(fun, lb, ub, options);
            end
            for i = n(1)+1:n(1)+n(2)
                % 用户i的本地优化问题
                z0=sum(x)-x(i)-(G+obj.PV(t)+B);
                fun = @(xi) user_local_optimization(xi, D(i), lambda, z0, ...
                    obj.beta, obj.sigma, rho_t, SOC_ev(i-n(1)));

                % 边界约束
                lb = 0;
                ub = 1 * D(i); % 最大不超过2.5倍基准需求

                % 求解一维优化问题
                options = optimset('Display', 'off', 'TolX', 1e-6);
                x(i) = fminbnd(fun, lb, ub, options);
            end
        end
         
        function [G,B,z,Tab,Tac,Tbc,conver_mu]=argmin_other2(obj,x1,x2,x3,lambda, rho_t,t,SOC_current) 
            % 系统优化：最小化电网成本+电池成本+ADMM惩罚项
            tmu=[0;0;0];conver_mu=zeros(2,20);
            for k2=1:40
                [G_new, B_new, Tab, Tac] = system_optimizationA(...
                    obj.PV(t), SOC_current(1), ...
                    lambda(1), sum(x1), rho_t(1), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);
                [G_newB, B_newB, bTab, Tbc] = system_optimizationB(...
                    obj.PV(t), SOC_current(2), ...
                    lambda(2), sum(x2), rho_t(2), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);
                [G_newC, B_newC, bTac, bTbc] = system_optimizationC(...
                    obj.PV(t), SOC_current(3), ...
                    lambda(3), sum(x3), rho_t(3), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);

                mu_old=tmu;
                primal_res=[Tab-bTab;Tac-bTac;Tbc-bTbc];
                tmu=tmu-0.1*primal_res;

                err1=norm(primal_res);
                err2=norm(tmu-mu_old); 
                conver_mu(:,k2)=[err1;err2]; 

                disp(['mu err:',num2str([err1,err2])])
                disp(['mu :',num2str(tmu')])
                disp(['3T:',num2str([Tab,Tac,Tbc])])

                if err1<1e-3 && err2 <1e-3
                    disp(k2)
                    
                    break;
                end
            end
            G = [G_new;G_newB;G_newC];
            B = [B_new;B_newB;B_newC];
            z=G+obj.PV(:,t)+B+[-Tab-Tac;bTab-Tbc;bTac+ bTbc]; 
        end

        function [u_utility,g_cost,b_cost,T_cost]=cost_compute(obj,x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current)
            u_utility = [0;0;0];  
            for i = 1:obj.n(1)
                u_utility(1) = u_utility(1) + piecewise_log_utility(...
                    x1(i), obj.D1(i, t), obj.beta, obj.sigma);
            end
            for i = 1:obj.n(2)
                u_utility(2) = u_utility(2) + piecewise_log_utility(...
                    x2(i), obj.D2(i, t), obj.beta, obj.sigma);
            end
            for i = 1:obj.n(3)
                u_utility(3) = u_utility(3) + piecewise_log_utility(...
                    x3(i), obj.D3(i, t), obj.beta, obj.sigma);
            end
            C_grid=@(x) obj.a0*x^2 + obj.b0*x + obj.c0;

            g_cost = [C_grid(G(1));C_grid(G(2));C_grid(G(3))];
            b_cost  = [battery_cost(B(1), SOC_current(1));battery_cost(B(2), SOC_current(2));battery_cost(B(3), SOC_current(3))];

            T_cost=obj.CofT*[Tab^2+Tac^2;Tab^2+Tbc^2;Tac^2+Tbc^2];
        end
        
        function draw(obj)
            close all;

            %% 可视化结果
            visualize_results(obj.X1_opt, obj.G_opt(1,:), obj.B_opt(1,:), obj.PV(1,:), obj.SOC_opt(1,:), obj.D1, ...
                obj.covnerge_numer, obj.total_cost_all(1),obj.Lambda_opt(1,:),-min(obj.Tab_opt,0),-min(obj.Tac_opt,0),'b2a','c2a',max(obj.Tab_opt,0),max(obj.Tac_opt,0),'a2b','a2c'); 
            visualize_results(obj.X2_opt, obj.G_opt(2,:), obj.B_opt(2,:), obj.PV(2,:), obj.SOC_opt(2,:), obj.D2, ...
                obj.covnerge_numer, obj.total_cost_all(2),obj.Lambda_opt(2,:),max(obj.Tab_opt,0),-min(obj.Tbc_opt,0),'a2b','c2b',-min(obj.Tab_opt,0),max(obj.Tbc_opt,0),'b2a','b2c');
            visualize_results(obj.X3_opt, obj.G_opt(3,:), obj.B_opt(3,:), obj.PV(3,:), obj.SOC_opt(3,:), obj.D3, ...
                obj.covnerge_numer, obj.total_cost_all(3),obj.Lambda_opt(3,:),max(obj.Tac_opt,0),max(obj.Tbc_opt,0),'a2c','b2c',-min(obj.Tac_opt,0),-min(obj.Tbc_opt,0),'c2a','c2b');


            %% 保存结果
            % save('optimization_results.mat', 'X_opt', 'G_opt', 'B_opt', 'SOC_opt', ...
            % 'PV', 'D','total_cost_all', 'covnerge_numer');


              % Display comprehensive power transfer situation of six subplots
            figure(4)  % Create large figure window

            % Define color scheme - using more professional colors
            colors = [0.2, 0.4, 0.8;    % Blue - Microgrid a
                0.8, 0.2, 0.2;    % Red - Microgrid b
                0.2, 0.7, 0.3];   % Green - Microgrid c

            lineStyles = {'-', '--', ':'};
            lineWidth = 2;

            % Subplot 1: Power transfer between three microgrids 
            x = 1:obj.T;

            h1 = plot(x, obj.Tab_opt, 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.Tac_opt, 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.Tbc_opt, 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Transfer Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(a) Inter-microgrid Power Transfer', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MGa2b', 'MGa2c', 'MG b2c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')


             figure(5)  % Create large figure window 
            x = 1:obj.T;

            plot(x, sum(obj.D1,1)-obj.PV(1,:), 'r-', 'LineWidth', 2, 'Marker', 'none');
            hold on
           plot(x, sum(obj.D2,1)-obj.PV(2,:), 'b-', 'LineWidth', 2, 'Marker', 'none');
           plot(x, sum(obj.D3,1)-obj.PV(3,:), 'g-', 'LineWidth', 2, 'Marker', 'none');

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Transfer Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(a) Inter-microgrid Power Transfer', 'FontSize', 12, 'FontWeight', 'bold')
            legend( 'MGA', 'MGB', 'MGC' )
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

        end

        
    end
    methods (Static)

    end
end
function cost = user_local_optimization(xi, Di, lambda, z, beta, sigma, rho,SOC)
% 用户本地优化目标函数
if ischar(SOC)
    utility = piecewise_log_utility(xi, Di, beta, sigma);
else
    utility = EV_utility(xi, Di, SOC); 
end
consensus_term = lambda * (xi) + (rho/2) * (xi + z)^2;
cost = -utility + consensus_term;
end
function U = piecewise_log_utility(x, D, beta, sigma)
% 分段对数效用函数
if x <= D
    U = beta * log(1 + sigma * x);
else
    U = beta * log(1 + sigma * D);
end
end
function U = EV_utility(x, D, SOC)
beta_ev=8*(1-SOC);
omega_ev=1;
% 分段对数效用函数
U=beta_ev*log(omega_ev*min(x,D)+1)/log(3);

end
function cost = battery_cost(x, SOC)
% 电池成本函数
% x: 电池功率 (kW), 正为放电，负为充电
% SOC: 当前荷电状态

if x <= 0  % 充电
    % (1-SOC)/SOC * exp(3) * 1 * (exp(x) - 1)
    cost = (1-SOC)/SOC * exp(0) * 0.01 * (exp(x) - 1);
else  % 放电
    % 0.01*x^2 + (1-SOC)/SOC * exp(3) * 1 * x
    cost = 0.00011*x^2 + (1-SOC)/SOC * exp(0) * 0.01 * x;
end
end 

function [SOC_new,SOC_ev_new] = update_battery_SOC(SOC_old, B, E_max, SOC_min, SOC_max,SOC_ev_old,Pev,E_ev)
% 更新电池SOC
% B: 电池功率 (kW), 正为放电，负为充电
% 假设时间间隔为1小时

eta_ch = 0.95;  % 充电效率
eta_dis = 0.97; % 放电效率

delta_SOC=zeros(3,1);
for kk=1:3
if B(kk) > 0  % 放电
    delta_SOC(kk) = -B(kk) / (eta_dis * E_max);
else  % 充电
    delta_SOC(kk) = -B(kk) * eta_ch / E_max;
end
end

SOC_new = SOC_old + delta_SOC;
SOC_new = max(SOC_min, min(SOC_max, SOC_new));

SOC_ev_new=SOC_ev_old+Pev/E_ev;
SOC_ev_new =  min(1, SOC_ev_new);
end


function visualize_results(X, G, B, PV, SOC, D, ...
    covnerge_numer, ~, lambda, input1, input2, input1T, input2T, output1, output2, output1T, output2T)
% Visualize all results


[~, T] = size(X);

figure('Position', [100, 100, 1400, 900]);

% 1. Power balance chart
subplot(3, 3, 1);
sum_x = sum(X, 1);
total_supply = G + PV + max(0, B); % Battery discharge as positive

plot(1:T, sum_x, 'b-', 'LineWidth', 2); hold on;
plot(1:T, total_supply, 'r--', 'LineWidth', 2);
plot(1:T, G, 'g:', 'LineWidth', 1.5);
plot(1:T, PV, 'y:', 'LineWidth', 1.5);

xlabel('Time (hour)');
ylabel('Power (kW)');
title('Power Balance');
legend('Total Demand', 'Total Supply', 'Grid', 'PV', 'Location', 'best');
grid on;
xlim([1, T]);

% 2. Power supply composition (stacked area chart)
subplot(3, 3, 2);
supply_components = [PV; G; max(0, B); input1; input2]';
area(1:T, supply_components);
xlabel('Time (hour)');
ylabel('Power (kW)');
title('Power Supply Composition');
legend({'PV', 'Grid', 'Battery Discharge', input1T, input2T}, 'Location', 'best');
grid on;
xlim([1, T]);

% 3. Power consumption composition
subplot(3, 3, 3);
supply_components = [output1; output2; -min(0, B)]';
area(1:T, supply_components);
xlabel('Time (hour)');
ylabel('Power (kW)');
title('Power Consumption Composition');
legend({output1T, output2T, 'Battery Charging'}, 'Location', 'best');
grid on;
xlim([1, T]);

% 4. Battery status
subplot(3, 3, 4);
yyaxis left;
plot(1:T, B, 'b-', 'LineWidth', 2);
ylabel('Charge/Discharge Power (kW)');
ylim([-max(abs(B))*1.1, max(abs(B))*1.1]);

yyaxis right;
plot(1:T, SOC, 'r-', 'LineWidth', 2);
ylabel('SOC');
ylim([0, 1]);

xlabel('Time (hour)');
title('Battery Status');
grid on;
xlim([1, T]);

% 5. User power consumption distribution
subplot(3, 3, 5);
imagesc(X);
colorbar;
xlabel('Time (hour)');
ylabel('User');
title('User Power Consumption Distribution (kW)');
colormap('hot');

% 6. Cost function curve
subplot(3, 3, 6);

% Plot battery cost function curve (example)
SOC_test = 0.3:0.15:0.9;
B_test = -20:2:20;
cost_matrix = zeros(length(SOC_test), length(B_test));

for i = 1:length(SOC_test)
    for j = 1:length(B_test)
        cost_matrix(i, j) = battery_cost(B_test(j), SOC_test(i));
    end
end

surf(B_test, SOC_test, cost_matrix);
xlabel('Battery Power (kW)');
ylabel('SOC');
zlabel('Cost (¥)');
title('Battery Cost Function');
grid on;

% 7. ADMM convergence
subplot(3, 3, 7);
bar(1:T, covnerge_numer(:, 1));
xlabel('Time Period');
ylabel('ADMM Iterations');
title('ADMM Convergence Speed by Period');
grid on;
ylim([0, max(covnerge_numer(:, 1))*1.1]);

% 8. Dual variable
subplot(3, 3, 8);
plot(1:T, lambda, 'b-', 'LineWidth', 2); hold on;
xlabel('Time (hour)');
ylabel('Value');
legend('Dual variable', 'Location', 'best');
grid on;
xlim([1, T]);

% 9. Baseline demand vs optimized demand
subplot(3, 3, 9);
baseline_demand = sum(D, 1);
optimized_demand = sum(X, 1);

plot(1:T, baseline_demand, 'b-', 'LineWidth', 2); hold on;
plot(1:T, optimized_demand, 'r--', 'LineWidth', 2);

xlabel('Time (hour)');
ylabel('Total Power Consumption (kW)');
title('Baseline Demand vs Optimized Demand');
legend('Baseline Demand', 'Optimized Demand', 'Location', 'best');
grid on;
xlim([1, T]);

% Add overall title
sgtitle('Electric Power System Optimization Dispatch Results Analysis', 'FontSize', 14, 'FontWeight', 'bold');
end

function [D, PV] = generationA(n, T)
 
A=getData(n);
PV=sum(A.GG,1); 
D=A.GC;
end

% function [D, PV] = generationA(n, T) 
% A=getData(n); 
% load('load_data.mat')
% D=resedential_load(25:48)/100;
% 
% PV=sum(A.GG,1); 
% PV=PV/max(PV)*2*max(D);
% end
function [D, PV] = generationB(n, T) 
% Generate user baseline demand 
load('load_data.mat')
D=commertial_load(25:48)/100;
PV = zeros(1, T);
daylight_hours = 6:19; % Daylight hours
sun_power = 25*3; % Peak PV power
PV(daylight_hours) = sun_power * sin(pi*(daylight_hours-6)/(19-6)).^2;
PV = PV + min(0, 2*randn(1, T)); % Add noise
PV = max(0, PV)/max(PV)*2*max(D); % Ensure non-negative
end 

function [D, PV] = generationC(n, T) 
% Generate user baseline demand
load('load_data.mat')
D=industrial_load(25:48)/10; 
PV = zeros(1, T);
daylight_hours = 6:19; % Daylight hours
sun_power = 25*3; % Peak PV power
PV(daylight_hours) = sun_power * sin(pi*(daylight_hours-6)/(19-6)).^2;
PV = PV + min(0, 2*randn(1, T)); % Add noise
%PV = max(0, PV)/max(PV)*2*max(D); % Ensure non-negative
PV = max(0, PV); % Ensure non-negative 
end