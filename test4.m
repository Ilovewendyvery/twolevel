%% 基于ADMM的电力系统优化调度程序
% 成本函数：
% C_grid(x) = a0*x^2 + b0*x + c0
% C_EBSS(x,SOC) = (1-SOC)/SOC*exp(3)*0.1*(exp(x)-1) .* (x <= 0) 
%                + (0.01*x^2 + (1-SOC)/SOC*exp(3)*0.1*x) .* (x > 0)

clear all; close all; clc;

%% 参数设置
n = 30;           % 用户数量
T = 24;          % 时间周期数（24小时）
beta = 1.5;      % 效用函数系数
sigma = 0.8;     % 效用函数形状参数

% 电网成本参数
a0 = 0.001;       % 二次项系数 (¥/kW²)
b0 = 0.02;        % 一次项系数 (¥/kW)  
c0 = 0;          % 常数项

% 电池参数
E_max = 1000;     % 电池容量 (kWh)
B_max = 250;      % 最大充放电功率 (kW)
SOC_min = 0.15;  % 最小SOC
SOC_max = 0.95;  % 最大SOC
SOC_initial = 0.5; % 初始SOC

% ADMM参数
rho = 0.8;       % 初始惩罚参数
max_admm_iter = 200; % 每个时段最大ADMM迭代次数
abs_tol = 1e-5;
rel_tol = 1e-5;

%% 生成模拟数据
rng(2023); % 设置随机种子

% 生成用户基准需求 D^A_{i,t}
D = zeros(n, T);
for i = 1:n
    % 每个用户有不同的用电模式
    user_type = randi([1, 3]); % 1:住宅, 2:商业, 3:工业
    base_load = [1.5, 3.0, 8.0]; % 基础负荷
    
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
sun_power = 25*7; % 峰值光伏功率
PV(daylight_hours) = sun_power * sin(pi*(daylight_hours-6)/(19-6)).^2;
PV = PV + 2*randn(1, T); % 添加噪声
PV = max(0, PV); % 确保非负

% 分时电价（用于电网成本中的b0调整）
time_price = zeros(1, T);
peak_hours = [8:11, 18:21];    % 峰时
flat_hours = [7, 12:17, 22];   % 平时
valley_hours = [1:6, 23:24];   % 谷时
time_price(peak_hours) = 0.8;  % 峰时电价
time_price(flat_hours) = 0.5;  % 平时电价
time_price(valley_hours) = 0.3; % 谷时电价

%% 存储结果
X_opt = zeros(n, T);      % 最优用电量 [kW]
G_opt = zeros(1, T);      % 最优电网供电 [kW]
B_opt = zeros(1, T);      % 最优电池充放电 [kW]（正为放电，负为充电）
SOC_opt = zeros(1, T);    % SOC轨迹
Lambda_opt = zeros(1, T); % 最优对偶变量

%% 成本函数定义
% 电网成本函数
C_grid_func = @(x) a0*x^2 + b0*x + c0;

% 电池成本函数（根据提供的公式）
C_battery_func = @(x, SOC) battery_cost(x, SOC);

% 总成本函数（用于评估）
total_cost_func = @(X, G, B, SOC) calculate_total_cost(X, G, B, SOC, D, ...
    beta, sigma, a0, b0, c0, time_price);

%% 初始化SOC
SOC_current = SOC_initial;

%% 主循环：按时间顺序处理（时间层在外层）
total_cost_all = 0;
convergence_info = zeros(T, 2); % 存储每个时段的收敛信息

for t = 1:T
    fprintf('\n=== 时段 %02d:00 - %02d:00 ===\n', t-1, t);
    fprintf('光伏预测: %.2f kW, 分时电价: %.2f ¥/kWh\n', PV(t), time_price(t));
    
    % 当前时段数据
    D_t = D(:, t);
    PV_t = PV(t);
    time_price_t = time_price(t);
    
    %% ADMM初始化（当前时段）
    rho_t = rho;  % 当前时段的惩罚参数
    x = D_t;      % 初始化为基准需求
    G = max(0, sum(x) - PV_t); % 初始电网供电
    B = 0;        % 初始电池功率 
    lambda = 0;   % 对偶变量
    z_prev=G+PV_t+B;
    
    %% 存储迭代历史（用于分析和调试）
    primal_res_history = zeros(max_admm_iter, 1);
    dual_res_history = zeros(max_admm_iter, 1);
    cost_history = zeros(max_admm_iter, 1);
    
    %% ADMM迭代（当前时段）
    converged = false;
    iter_used = max_admm_iter;
    
    for k = 1:max_admm_iter
        % 保存上一次的z值（用于计算对偶残差）
        if k > 1
            z_prev = z;
        end
        
        %% 步骤1: 更新x（用户侧，可并行）
        total_x_prev = sum(x);
        
        % 并行更新每个用户的用电量
        for i = 1:n  
            % 用户i的本地优化问题
            z0=sum(x)-x(i)-(G+PV_t+B);
            fun = @(xi) user_local_optimization(xi, D_t(i), lambda, z0, ...
                beta, sigma, rho_t);
            
            % 边界约束
            lb = 0;
            ub = 1 * D_t(i); % 最大不超过2.5倍基准需求
            
            % 求解一维优化问题
            options = optimset('Display', 'off', 'TolX', 1e-6);
            x(i) = fminbnd(fun, lb, ub, options);
        end
        
        %% 步骤2: 更新系统变量 (z, G, B) 
        
        % 系统优化：最小化电网成本+电池成本+ADMM惩罚项
        [G_new, B_new, feasible, system_cost_val] = system_optimization(...
             PV_t, SOC_current, ...
            lambda, sum(x), rho_t, time_price_t, ...
            a0, b0, B_max, E_max);
        
        if feasible
            G = G_new;
            B = B_new; 
            z=G+PV_t+B;  
        else
            warning('时段 %d 迭代 %d: 系统问题不可行，保持上一解', t, k);
        end
        
        %% 步骤3: 更新对偶变量
        primal_res = sum(x) - (G+PV_t+B);
        lambda = lambda + rho_t * primal_res;
        
        %% 计算收敛指标
        primal_res_history(k) = abs(primal_res);
        
        if k > 1
            dual_res_history(k) = rho_t * abs(z - z_prev);
        else
            dual_res_history(k) = Inf;
        end
        
        % 计算当前解的总成本（用于监控）
        user_utility = 0;
        for i = 1:n
            user_utility = user_utility + piecewise_log_utility(...
                x(i), D_t(i), beta, sigma);
        end
        
        grid_cost_val = C_grid_func(G) + time_price_t * G;
        battery_cost_val = C_battery_func(B, SOC_current); 
        
        cost_history(k) = -user_utility + grid_cost_val + ...
            battery_cost_val;
        
        %% 收敛检查
        eps_primal = abs_tol + rel_tol * max(abs(sum(x)), abs(G+PV_t+B));
        eps_dual = abs_tol + rel_tol * abs(z);
        
        if k > 1 && primal_res_history(k) < eps_primal && ...
                dual_res_history(k) < eps_dual
            fprintf('  ADMM收敛于 %d 次迭代\n', k);
            fprintf('  原始残差: %.2e, 对偶残差: %.2e\n', ...
                primal_res_history(k), dual_res_history(k));
            converged = true;
            iter_used = k;
            break;
        end
        
        %% 自适应调整rho（每10次迭代）
        if mod(k, 10) == 0 && k > 10
            mu = 10;  % 调整阈值
            tau = 1.5; % 调整因子
            
            if primal_res_history(k) > mu * dual_res_history(k)
                rho_t = rho_t * tau;
                fprintf('  迭代 %d: 增加rho到 %.3f\n', k, rho_t);
            elseif dual_res_history(k) > mu * primal_res_history(k)
                rho_t = rho_t / tau;
                fprintf('  迭代 %d: 减少rho到 %.3f\n', k, rho_t);
            end
        end
        
        %% 显示进度
        if mod(k, 20) == 0 && k < max_admm_iter
            fprintf('  迭代 %3d: 原始残差=%.2e, 对偶残差=%.2e, 成本=%.2f\n', ...
                k, primal_res_history(k), dual_res_history(k), cost_history(k));
        end
    end
    
    if ~converged
        fprintf('  ADMM未收敛，使用最后迭代结果\n');
    end
    
    %% 存储当前时段最优解
    X_opt(:, t) = x;
    G_opt(t) = G;
    B_opt(t) = B;
    SOC_opt(t) = SOC_current;
    Lambda_opt(t) = lambda;
    
    % 记录收敛信息
    convergence_info(t, 1) = iter_used;
    convergence_info(t, 2) = primal_res_history(min(iter_used, max_admm_iter));
    
    % 计算当前时段成本
    user_utility_t = 0;
    for i = 1:n
        user_utility_t = user_utility_t + piecewise_log_utility(...
            x(i), D_t(i), beta, sigma);
    end
    
    grid_cost_t = C_grid_func(G) + time_price_t * G;
    battery_cost_t = C_battery_func(B, SOC_current);
    total_cost_t = -user_utility_t + grid_cost_t + battery_cost_t;
    total_cost_all = total_cost_all + total_cost_t;
    
    fprintf('  优化结果: 总需求=%.2f kW, 电网=%.2f kW, 电池=%.2f kW, SOC=%.3f\n', ...
        sum(x), G, B, SOC_current);
    fprintf('  时段成本: 用户效用=%.2f, 电网成本=%.2f, 电池成本=%.2f, 总成本=%.2f\n', ...
        user_utility_t, grid_cost_t, battery_cost_t, total_cost_t);

    % 更新SOC（供下一个时段使用）
            SOC_current = update_battery_SOC(SOC_current, B_new, ...
                E_max, SOC_min, SOC_max);
    
end

%% 计算总体性能指标
fprintf('\n\n=== 优化完成 ===\n');
fprintf('总成本: %.2f ¥\n', total_cost_all);
fprintf('平均每时段成本: %.2f ¥\n', total_cost_all/T);
fprintf('总用电量: %.2f kWh\n', sum(sum(X_opt)));
fprintf('总电网供电: %.2f kWh\n', sum(G_opt));
fprintf('总电池充电: %.2f kWh\n', sum(B_opt(B_opt < 0)));
fprintf('总电池放电: %.2f kWh\n', sum(B_opt(B_opt > 0)));
fprintf('最终SOC: %.3f\n', SOC_current);
fprintf('平均ADMM迭代次数: %.1f\n', mean(convergence_info(:,1)));

%% 成本分解分析
[total_utility, total_grid_cost, total_battery_cost] = ...
    analyze_costs(X_opt, G_opt, B_opt, SOC_opt, D, ...
    beta, sigma, a0, b0, c0, time_price);

fprintf('\n=== 成本分解 ===\n');
fprintf('用户效用总价值: %.2f ¥\n', total_utility);
fprintf('电网总成本: %.2f ¥ (%.1f%%)\n', ...
    total_grid_cost, 100*total_grid_cost/total_cost_all);
fprintf('电池总成本: %.2f ¥ (%.1f%%)\n', ...
    total_battery_cost, 100*total_battery_cost/total_cost_all);
fprintf('净总成本: %.2f ¥\n', -total_utility + total_grid_cost + total_battery_cost);

%% 可视化结果
visualize_results(X_opt, G_opt, B_opt, PV, SOC_opt, D, time_price, ...
    convergence_info, total_cost_all,Lambda_opt);

%% 保存结果
save('optimization_results.mat', 'X_opt', 'G_opt', 'B_opt', 'SOC_opt', ...
    'PV', 'D', 'time_price', 'total_cost_all', 'convergence_info');

%% ==================== 辅助函数 ====================
function cost = battery_cost(x, SOC)
    % 电池成本函数
    % x: 电池功率 (kW), 正为放电，负为充电
    % SOC: 当前荷电状态
    
    if x <= 0  % 充电
        % (1-SOC)/SOC * exp(3) * 1 * (exp(x) - 1)
        cost = (1-SOC)/SOC * exp(0) * 0.1 * (exp(x) - 1);
    else  % 放电
        % 0.01*x^2 + (1-SOC)/SOC * exp(3) * 1 * x
        cost = 0.00011*x^2 + (1-SOC)/SOC * exp(0) * 0.1 * x;
    end
end

function U = piecewise_log_utility(x, D, beta, sigma)
    % 分段对数效用函数
    if x <= D
        U = beta * log(1 + sigma * x);
    else
        U = beta * log(1 + sigma * D);
    end
end

function cost = user_local_optimization(xi, Di, lambda, z, beta, sigma, rho)
    % 用户本地优化目标函数
    utility = piecewise_log_utility(xi, Di, beta, sigma);
    consensus_term = lambda * (xi) + (rho/2) * (xi + z)^2;
    cost = -utility + consensus_term;
end

function [G_opt, B_opt, feasible, total_cost] = system_optimization(...
     PV, SOC, lambda, sum_x, rho, time_price, ...
    a0, b0, B_max, E_max)
    % 系统优化：求解最优的G和B
    
    % 净需求
    net_demand = sum_x - PV;
    
    % 定义优化问题
    options = optimoptions('fmincon', 'Display', 'off', ...
        'Algorithm', 'sqp', 'MaxIterations', 100);
    
    % 初始猜测
    x0 = [max(0, net_demand), 0];  % [G, B]
    
    % 边界约束
    lb = [0, -B_max];
    ub = [Inf, B_max];
    
    % 线性等式约束：G + B = net_demand
    Aeq = [1, 1];
    beq = net_demand;
    
    % 目标函数
    function f = system_objective(vars)
        G = vars(1);
        B = vars(2);
        
        % 电网成本（考虑分时电价）
        grid_cost = a0*G^2 + (b0 + 0.05*time_price)*G;
        
        % 电池成本
        battery_cost_val = battery_cost(B, SOC);
        
        % ADMM惩罚项  
        admm_penalty = lambda * (-G-PV-B) + ...
            (rho/2) * (sum_x - (G+PV+B))^2;
        
        f = grid_cost + battery_cost_val + admm_penalty;
    end
    
    % 求解
    try
        [vars_opt, fval] = fmincon(@system_objective, x0, [], [], ...
            [], [], lb, ub, [], options);
        
        G_opt = vars_opt(1);
        B_opt = vars_opt(2);
        total_cost = fval;
        feasible = true;
    catch
        % 如果优化失败，使用简单启发式
        G_opt = max(0, net_demand);
        B_opt = 0;
        total_cost = Inf;
        feasible = false;
    end
end

function SOC_new = update_battery_SOC(SOC_old, B, E_max, SOC_min, SOC_max)
    % 更新电池SOC
    % B: 电池功率 (kW), 正为放电，负为充电
    % 假设时间间隔为1小时
    
    eta_ch = 0.95;  % 充电效率
    eta_dis = 0.97; % 放电效率
    
    if B > 0  % 放电
        delta_SOC = -B / (eta_dis * E_max);
    else  % 充电
        delta_SOC = -B * eta_ch / E_max;
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

function total_cost = calculate_total_cost(X, G, B, SOC, D, ...
    beta, sigma, a0, b0, c0, time_price)
    % 计算总成本
    
    [total_utility, total_grid_cost, total_battery_cost] = ...
        analyze_costs(X, G, B, SOC, D, beta, sigma, a0, b0, c0, time_price);
    
    total_cost = -total_utility + total_grid_cost + total_battery_cost;
end

function visualize_results(X, G, B, PV, SOC, D, time_price, ...
    convergence_info, total_cost,lambda)
    % 可视化所有结果
    
    [n, T] = size(X);
    
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
    bar(1:T, convergence_info(:,1));
    xlabel('时段');
    ylabel('ADMM迭代次数');
    title('各时段ADMM收敛速度');
    grid on;
    ylim([0, max(convergence_info(:,1))*1.1]);
    
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