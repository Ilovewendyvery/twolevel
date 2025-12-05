classdef New < handle
    properties
        %% 参数设置
        n;           % 每个微网用户数量
        T = 24;          % 时间周期数（24小时）
        beta = 15;      % 效用函数系数
        sigma = 8;     % 效用函数形状参数

        % 电网成本参数
        a0 = 0.1;       % 二次项系数 (¥/kW²)
        b0 = 2;        % 一次项系数 (¥/kW)
        c0 = 0;          % 常数项

        % 电池参数
        E_max = 1000;     % 电池容量 (kWh)
        B_max = 20;      % 最大充放电功率 (kW)
        SOC_min = 0.15;  % 最小SOC
        SOC_max = 0.95;  % 最大SOC
        SOC_initial = 0.5; % 初始SOC

        %传递参数
        CofT=0.1;
        Tmax=10;

        % ADMM参数
        rho = 0.8;       % 初始惩罚参数
        max_admm_iter = 200; % 每个时段最大ADMM迭代次数
        abs_tol = 1e-5;
        rel_tol = 1e-5;


        PV;        %(3,T)
        time_price;%(1,T)
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
        converge_primal_lam_end;%(T,3)
        converge_primal_mu_end;%(T,3)
    end

    methods
        function obj = New()
            %生成模拟数据
            obj.T=24; obj.n=[30,30,30];
            rng(2023); % 设置随机种子 

            %% For A B C
            [D,PV]=generation(obj.n(1),obj.T,1);obj.D1=D;obj.PV(1,:)=PV;
            [D,PV]=generation(obj.n(2),obj.T,2);obj.D2=D;obj.PV(2,:)=PV;
            [D,PV]=generation(obj.n(3),obj.T,3);obj.D3=D;obj.PV(3,:)=PV;  

            % 分时电价（用于电网成本中的b0调整）
            time_price = zeros(1, obj.T);
            peak_hours = [8:11, 18:21];    % 峰时
            flat_hours = [7, 12:17, 22];   % 平时
            valley_hours = [1:6, 23:24];   % 谷时
            time_price(peak_hours) = 0.8;  % 峰时电价
            time_price(flat_hours) = 0.5;  % 平时电价
            time_price(valley_hours) = 0.3; % 谷时电价
            obj.time_price=time_price;


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
            obj.converge_primal_lam_end= zeros(obj.T, 3);%(T,3)
            obj.converge_primal_mu_end=zeros(obj.T, 3);%(T,3) 
            obj.total_cost_all = 0;
        end

        function run(obj)
            tic;
            run_inital(obj); 
            SOC_current = obj.SOC_initial*ones(3,1);
            %% 主循环：按时间顺序处理（时间层在外层）
            for t = 1:obj.T
                obj.hour(t,SOC_current);

                % 更新SOC（供下一个时段使用）
                SOC_current = update_battery_SOC(SOC_current, obj.B_opt(:,t), ...
                    obj.E_max, obj.SOC_min, obj.SOC_max);
            end 
            toc; 
        end
        
        function hour(obj,t,SOC_current)
            fprintf('\n=== 时段 %02d:00 - %02d:00 ===\n', t-1, t);
            fprintf('光伏预测: %.2f kW, 分时电价: %.2f ¥/kWh\n', obj.PV(1,t), obj.time_price(t));

            % 当前时段数据 
            time_price_t = obj.time_price(t);

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
                x1=argmin_x(obj,x1,G(1),B(1),lambda(1), rho_t(1),t,obj.n(1),obj.D1(:, t));
                x2=argmin_x(obj,x2,G(2),B(2),lambda(2), rho_t(2),t,obj.n(2),obj.D2(:, t));
                x3=argmin_x(obj,x3,G(3),B(3),lambda(3), rho_t(3),t,obj.n(3),obj.D3(:, t));

                 

                %% 步骤2: 更新系统变量 (z, G, B)
                [G,B,z,Tab,Tac,Tbc]=argmin_other2(obj,x1,x2,x3,lambda, rho_t,t,SOC_current);



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
                [u_utility,g_cost,b_cost]=obj.cost_compute(x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current,time_price_t);
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
            % obj.covnerge_numer(t, 2) = max(primal_res_history(:,min(iter_used, obj.max_admm_iter)));

            % 计算当前时段成本

            [u_utility,g_cost,b_cost]=obj.cost_compute(x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current,time_price_t);
            cost_history(:,k) = -u_utility + g_cost + b_cost;

            obj.total_cost_all = obj.total_cost_all + cost_history(:,k);

            fprintf('  优化结果: 总需求=%.2f kW, 电网=%.2f kW, 电池=%.2f kW, SOC=%.3f\n', ...
                sum(x1), G(1), B(1), SOC_current(1));
            fprintf('  时段成本: 用户效用=%.2f, 电网成本=%.2f, 电池成本=%.2f, 总成本=%.2f\n', ...
                u_utility(1), g_cost(1), b_cost(1), obj.total_cost_all(1)); 

        end

        function x=argmin_x(obj,x,G,B,lambda, rho_t,t,n,D)
            % 并行更新每个用户的用电量
            for i = 1:n
                % 用户i的本地优化问题
                z0=sum(x)-x(i)-(G+obj.PV(t)+B);
                fun = @(xi) user_local_optimization(xi, D(i), lambda, z0, ...
                    obj.beta, obj.sigma, rho_t);

                % 边界约束
                lb = 0;
                ub = 1 * D(i); % 最大不超过2.5倍基准需求

                % 求解一维优化问题
                options = optimset('Display', 'off', 'TolX', 1e-6);
                x(i) = fminbnd(fun, lb, ub, options);
            end
        end
         
        function [G,B,z,Tab,Tac,Tbc]=argmin_other2(obj,x1,x2,x3,lambda, rho_t,t,SOC_current) 
            % 系统优化：最小化电网成本+电池成本+ADMM惩罚项
            tmu=[0;0;0];
            for k2=1:40
                [G_new, B_new, Tab, Tac] = system_optimizationA(...
                    obj.PV(t), SOC_current(1), ...
                    lambda(1), sum(x1), rho_t(1), obj.time_price(t), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);
                [G_newB, B_newB, bTab, Tbc] = system_optimizationB(...
                    obj.PV(t), SOC_current(2), ...
                    lambda(2), sum(x2), rho_t(2), obj.time_price(t), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);
                [G_newC, B_newC, bTac, bTbc] = system_optimizationC(...
                    obj.PV(t), SOC_current(3), ...
                    lambda(3), sum(x3), rho_t(3), obj.time_price(t), ...
                    obj.a0, obj.b0, obj.B_max, obj.E_max,tmu,obj.CofT,obj.Tmax);

                mu_old=tmu;
                primal_res=[Tab-bTab;Tac-bTac;Tbc-bTbc];
                tmu=tmu-0.1*primal_res;

                err1=norm(primal_res);
                err2=norm(tmu-mu_old);
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

        function [u_utility,g_cost,b_cost,T_cost]=cost_compute(obj,x1,x2,x3,t,G,B,Tab,Tac,Tbc,SOC_current,time_price_t)
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

            g_cost = [C_grid(G(1)) + time_price_t * G(1);C_grid(G(2)) + time_price_t * G(3);C_grid(G(3)) + time_price_t * G(3)];
            b_cost  = [battery_cost(B(1), SOC_current(1));battery_cost(B(2), SOC_current(2));battery_cost(B(3), SOC_current(3))];

            T_cost=obj.CofT*[Tab^2+Tac^2;Tab^2+Tbc^2;Tac^2+Tbc^2];
        end
        
        function draw(obj)
            %  %% 计算总体性能指标
            % fprintf('\n\n=== 优化完成 ===\n');
            % fprintf('总成本: %.2f ¥\n', obj.total_cost_all);
            % fprintf('平均每时段成本: %.2f ¥\n', obj.total_cost_all/obj.T);
            % fprintf('总用电量: %.2f kWh\n', sum(sum(obj.X1_opt)));
            % fprintf('总电网供电: %.2f kWh\n', sum(obj.G_opt(1,:)));
            % fprintf('总电池充电: %.2f kWh\n', sum(obj.B_opt(1,obj.B_opt(1,:) < 0)));
            % fprintf('总电池放电: %.2f kWh\n', sum(obj.B_opt(1,obj.B_opt(1,:) > 0)));
            % fprintf('最终SOC: %.3f\n', SOC_current);
            % fprintf('平均ADMM迭代次数: %.1f\n', mean(obj.covnerge_numer(:,1)));
            % 
            % %% 成本分解分析
            % [total_utility, total_grid_cost, total_battery_cost] = ...
            %     analyze_costs(obj.X1_opt, obj.G_opt, obj.B_opt, obj.SOC_opt, obj.D1, ...
            %     obj.beta, obj.sigma, obj.a0, obj.b0, obj.c0, obj.time_price);
            % 
            % fprintf('\n=== 成本分解 ===\n');
            % fprintf('用户效用总价值: %.2f ¥\n', total_utility);
            % fprintf('电网总成本: %.2f ¥ (%.1f%%)\n', ...
            %     total_grid_cost, 100*total_grid_cost/obj.total_cost_all);
            % fprintf('电池总成本: %.2f ¥ (%.1f%%)\n', ...
            %     total_battery_cost, 100*total_battery_cost/obj.total_cost_all);
            % fprintf('净总成本: %.2f ¥\n', -total_utility + total_grid_cost + total_battery_cost);
            % 

            close all;

            %% 可视化结果
            visualize_results(obj.X1_opt, obj.G_opt(1,:), obj.B_opt(1,:), obj.PV(1,:), obj.SOC_opt(1,:), obj.D1, obj.time_price, ...
                obj.covnerge_numer, obj.total_cost_all(1),obj.Lambda_opt(1,:));

            visualize_results(obj.X2_opt, obj.G_opt(2,:), obj.B_opt(2,:), obj.PV(2,:), obj.SOC_opt(2,:), obj.D2, obj.time_price, ...
                obj.covnerge_numer, obj.total_cost_all(2),obj.Lambda_opt(2,:));
            visualize_results(obj.X3_opt, obj.G_opt(3,:), obj.B_opt(3,:), obj.PV(3,:), obj.SOC_opt(3,:), obj.D3, obj.time_price, ...
                obj.covnerge_numer, obj.total_cost_all(3),obj.Lambda_opt(3,:));


            %% 保存结果
            % save('optimization_results.mat', 'X_opt', 'G_opt', 'B_opt', 'SOC_opt', ...
            % 'PV', 'D', 'time_price', 'total_cost_all', 'covnerge_numer');


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

        end

        
    end
    methods (Static)

    end
end
function cost = user_local_optimization(xi, Di, lambda, z, beta, sigma, rho)
% 用户本地优化目标函数
utility = piecewise_log_utility(xi, Di, beta, sigma);
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

function SOC_new = update_battery_SOC(SOC_old, B, E_max, SOC_min, SOC_max)
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
end

function [total_utility, total_grid_cost, total_battery_cost] = ...
    analyze_costs(X, G, B, SOC, D, beta, sigma, a0, b0, c0, time_price)
% 分析总成本分解

[n, T] = size(X);

total_utility = 0;
total_grid_cost = 0;
total_battery_cost = 0;

for t = 1:T
    % 用户效用
    for i = 1:n
        total_utility = total_utility + ...
            piecewise_log_utility(X(i,t), D(i,t), beta, sigma);
    end

    % 电网成本
    total_grid_cost = total_grid_cost + ...
        a0*G(t)^2 + (b0 + time_price(t))*G(t) + c0;

    % 电池成本
    total_battery_cost = total_battery_cost + ...
        battery_cost(B(t), SOC(t));
end
end

function visualize_results(X, G, B, PV, SOC, D, time_price, ...
    covnerge_numer, ~,lambda)
% 可视化所有结果

[~, T] = size(X);

figure('Position', [100, 100, 1400, 900]);

% 1. 功率平衡图
subplot(3, 3, 1);
sum_x = sum(X, 1);
total_supply = G + PV + max(0, B); % 电池放电为正

plot(1:T, sum_x, 'b-', 'LineWidth', 2); hold on;
plot(1:T, total_supply, 'r--', 'LineWidth', 2);
plot(1:T, G, 'g:', 'LineWidth', 1.5);
plot(1:T, PV, 'y:', 'LineWidth', 1.5);

xlabel('时间 (小时)');
ylabel('功率 (kW)');
title('功率平衡');
legend('总需求', '总供应', '电网', '光伏', 'Location', 'best');
grid on;
xlim([1, T]);

% 2. 电源组成（堆叠面积图）
subplot(3, 3, 2);
supply_components = [PV; G; max(0, B); -min(0, B)]';
area(1:T, supply_components);
xlabel('时间 (小时)');
ylabel('功率 (kW)');
title('电源组成');
legend({'光伏', '电网', '电池放电', '电池充电'}, 'Location', 'best');
grid on;
xlim([1, T]);

% 3. 电池状态
subplot(3, 3, 3);
yyaxis left;
plot(1:T, B, 'b-', 'LineWidth', 2);
ylabel('充放电功率 (kW)');
ylim([-max(abs(B))*1.1, max(abs(B))*1.1]);

yyaxis right;
plot(1:T, SOC, 'r-', 'LineWidth', 2);
ylabel('SOC');
ylim([0, 1]);

xlabel('时间 (小时)');
title('电池状态');
grid on;
xlim([1, T]);

% 4. 电价和光伏
subplot(3, 3, 4);
yyaxis left;
plot(1:T, time_price, 'b-', 'LineWidth', 2);
ylabel('电价 (¥/kWh)');

yyaxis right;
plot(1:T, PV, 'r-', 'LineWidth', 2);
ylabel('光伏功率 (kW)');

xlabel('时间 (小时)');
title('电价和光伏曲线');
grid on;
xlim([1, T]);

% 5. 用户用电量分布
subplot(3, 3, 5);
imagesc(X);
colorbar;
xlabel('时间 (小时)');
ylabel('用户');
title('用户用电量分布 (kW)');
colormap('hot');

% 6. 成本函数曲线
subplot(3, 3, 6);

% 绘制电池成本函数曲线（示例）
SOC_test = 0.3:0.1:0.9;
B_test = -20:0.5:20;
cost_matrix = zeros(length(SOC_test), length(B_test));

for i = 1:length(SOC_test)
    for j = 1:length(B_test)
        cost_matrix(i, j) = battery_cost(B_test(j), SOC_test(i));
    end
end

surf(B_test, SOC_test, cost_matrix);
xlabel('电池功率 (kW)');
ylabel('SOC');
zlabel('成本 (¥)');
title('电池成本函数');
grid on;

% 7. ADMM收敛情况
subplot(3, 3, 7);
bar(1:T, covnerge_numer(:,1));
xlabel('时段');
ylabel('ADMM迭代次数');
title('各时段ADMM收敛速度');
grid on;
ylim([0, max(covnerge_numer(:,1))*1.1]);

% % 8. 成本分解饼图
% subplot(3, 3, 8);
% [utility_val, grid_cost_val, battery_cost_val] = ...
%     analyze_costs(X, G, B, SOC, D, beta, sigma, a0, b0, c0, time_price);
%
% cost_components = [-utility_val, grid_cost_val, battery_cost_val];
% pie(cost_components);
% title(sprintf('总成本分解 (总计: %.1f¥)', total_cost));
% legend({'用户效用成本', '电网成本', '电池成本'}, 'Location', 'best');
subplot(3,3,8)
plot(1:T, lambda, 'b-', 'LineWidth', 2); hold on;

xlabel('时间 (小时)');
legend('price', 'Location', 'best');
grid on;
xlim([1, T]);


% 9. 基准需求 vs 优化需求
subplot(3, 3, 9);
baseline_demand = sum(D, 1);
optimized_demand = sum(X, 1);

plot(1:T, baseline_demand, 'b-', 'LineWidth', 2); hold on;
plot(1:T, optimized_demand, 'r--', 'LineWidth', 2);

xlabel('时间 (小时)');
ylabel('总用电量 (kW)');
title('基准需求 vs 优化需求');
legend('基准需求', '优化后需求', 'Location', 'best');
grid on;
xlim([1, T]);

% 添加总标题
sgtitle('电力系统优化调度结果分析', 'FontSize', 14, 'FontWeight', 'bold');
end

function [D,PV]=generation(n,T,user_type)
% 生成用户基准需求
D = zeros(n, T);
for i = 1:n
    % 每个用户有不同的用电模式
    %user_type = randi([1, 3]); % 1:住宅, 2:商业, 3:工业
    base_load = [2, 2, 2]; % 基础负荷
    %base_load = [1.5, 3.0, 8.0]; % 基础负荷


    % 生成日负荷曲线
    t_vec = 1:T;
    if user_type == 1  % 住宅用户：早晚高峰
        pattern = 0.7 + 0.3*sin(2*pi*(t_vec-7)/24) + ...
            0.4*exp(-((t_vec-8).^2)/(2*2^2)) + ...
            0.5*exp(-((t_vec-19).^2)/(2*3^2));
    elseif user_type == 2  % 商业用户：白天高峰
        pattern = 0.5 + 0.5*sin(2*pi*(t_vec-12)/24) + ...
            0.6*exp(-((t_vec-14).^2)/(2*4^2));
    else  % 工业用户：相对平稳
        pattern = 0.8 + 0.2*sin(2*pi*(t_vec)/24);
    end

    D(i, :) = base_load(user_type) * pattern + 0.2*randn(1, T);
    D(i, :) = max(0.1, D(i, :)); % 确保非负
end

% 光伏发电预测 PV_A
PV = zeros(1, T);
daylight_hours = 6:19; % 白天时段
sun_power = 25*3; % 峰值光伏功率
PV(daylight_hours) = sun_power * sin(pi*(daylight_hours-6)/(19-6)).^2;
PV = PV + 2*randn(1, T); % 添加噪声
PV = max(0, PV); % 确保非负
end