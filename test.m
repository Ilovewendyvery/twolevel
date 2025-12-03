%% 基于ADMM的电力系统优化调度程序（完整版）
clear all; close all; clc;

%% 参数设置
n = 20;           % 用户数量
T = 24;          % 时间周期数
beta = 1.0;      % 效用函数系数
sigma = 0.5;     % 效用函数形状参数

% 电网成本参数
grid_params.a = 0.01;    % 二次项系数 (¥/kW²)
grid_params.b = 0.5;     % 一次项系数 (¥/kW)
grid_params.c = 0;       % 常数项

% 电池成本参数
battery_params.cycle_cost_coeff = 0.02;   % 循环成本系数 (¥/kWh)
battery_params.SOC_penalty = 10;         % SOC偏离惩罚系数
battery_params.charge_loss = 0.05;       % 充电损失系数
battery_params.discharge_loss = 0.05;    % 放电损失系数
battery_params.B_max = 20;               % 最大充放电功率 (kW)
battery_params.SOC_min = 0.2;
battery_params.SOC_max = 0.9;
battery_params.SOC_target = 0.5;

% ADMM参数
rho = 1.0;       % 惩罚参数
max_iter = 1000;
abs_tol = 1e-4;
rel_tol = 1e-4;

%% 生成模拟数据
rng(42);
D = zeros(n, T);
for i = 1:n
    base_demand = 2 + 0.5*randn();
    daily_pattern = sin(2*pi*(0:T-1)/24) + 1;
    noise = 0.1*randn(1, T);
    D(i, :) = max(0.1, base_demand * daily_pattern + noise);
end

% 光伏发电
PV = zeros(1, T);
daylight_hours = 6:19;
PV(daylight_hours) = 15 * sin(pi*(daylight_hours-6)/(19-6));
PV = max(0, PV + 0.5*randn(1, T));

% 分时电价（模拟）
time_of_use_price = ones(1, T) * 0.5;  % 基础电价
time_of_use_price(7:10) = 0.8;        % 早高峰
time_of_use_price(18:21) = 0.8;       % 晚高峰
time_of_use_price(23:6) = 0.3;        % 谷时

%% 成本函数定义
% 电网成本函数（二次成本 + 分时电价）
C_grid_func = @(G, t) grid_params.a * G^2 + time_of_use_price(t) * G;

% 电池成本函数
C_battery_func = @(B, SOC, params) battery_total_cost(B, SOC, params);

%% 初始化变量
x = zeros(n, T);
G = zeros(1, T);
B = zeros(1, T);
lambda = zeros(1, T);
z = zeros(1, T);
x_tilde = zeros(n, T);

% SOC初始化
SOC = zeros(1, T);
SOC(1) = battery_params.SOC_target;



for t = 1:T
    %% ADMM主循环
    primal_residual = zeros(max_iter, 1);
    dual_residual = zeros(max_iter, 1);
    cost_history = zeros(max_iter, 1);
    for k = 1:max_iter
        %% 步骤1: 更新x（用户用电量）
        for i = 1:n

            % 用户i在时段t的子问题
            fun = @(xi) user_subproblem(xi, D(i,t), lambda(t), z(t), ...
                beta, sigma, rho, n);

            % 搜索范围 [0, 3*D]
            x_opt = fminbnd(fun, 0, 3*D(i,t));
            x_tilde(i, t) = x_opt;
        end

        %% 步骤2: 更新系统变量 (z, G, B)
        total_demand = sum(x_tilde(:, t));
        net_demand = total_demand - PV(t);

        % 如果是第一个时段，使用简单的优化
        if t == 1
            % 使用fmincon求解G和B
            options = optimoptions('fmincon', 'Display', 'off');

            % 初始猜测
            x0 = [net_demand/2, net_demand/2];

            % 界限约束
            lb = [0, -battery_params.B_max];
            ub = [Inf, battery_params.B_max];

            % 非线性约束：功率平衡
            nonlcon = @(vars) power_balance_constraint(vars, net_demand);

            % 目标函数
            obj_fun = @(vars) system_cost(vars(1), vars(2), ...
                SOC(t), t, C_grid_func, C_battery_func, battery_params);

            % 求解
            vars = fmincon(obj_fun, x0, [], [], [], [], lb, ub, ...
                nonlcon, options);

            G(t) = vars(1);
            B(t) = vars(2);
        else
            % 后续时段考虑SOC连续性
            SOC_prev = SOC(t-1);

            % 更新当前SOC（简化模型）
            delta_SOC = -B(t-1) / (battery_params.B_max * 4); % 假设容量为80kWh
            SOC_current = SOC_prev + delta_SOC;
            SOC_current = max(battery_params.SOC_min, ...
                min(battery_params.SOC_max, SOC_current));

            % 优化当前时段
            options = optimoptions('fmincon', 'Display', 'off');
            x0 = [net_demand/2, net_demand/2];

            % 根据SOC调整B的界限
            B_available = min(battery_params.B_max, ...
                (SOC_current - battery_params.SOC_min) * battery_params.B_max * 4);
            B_discharge_max = min(battery_params.B_max, B_available);

            B_charge_max = min(battery_params.B_max, ...
                (battery_params.SOC_max - SOC_current) * battery_params.B_max * 4);

            lb = [0, -B_charge_max];
            ub = [Inf, B_discharge_max];

            nonlcon = @(vars) power_balance_constraint(vars, net_demand);
            obj_fun = @(vars) system_cost(vars(1), vars(2), ...
                SOC_current, t, C_grid_func, C_battery_func, battery_params);

            vars = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, options);

            G(t) = vars(1);
            B(t) = vars(2);
            SOC(t) = SOC_current;
        end

        z(t) = total_demand; % 一致性变量

        %% 步骤3: 更新对偶变量
        primal_res = sum(x_tilde(:, t)) - z(t);
        lambda(t) = lambda(t) + rho * primal_res;

        %% 计算残差和成本
        primal_residual(k) = norm(sum(x_tilde(:, t), 1) - z(:, t), 2);

        if k > 1
            dual_residual(k) = rho * norm(z(:, t) - z_prev(:, t), 2);
        else
            dual_residual(k) = Inf;
        end

        % 计算总成本
        total_cost = 0;
        % 用户效用（负成本）
        utility_sum = 0;
        for i = 1:n
            utility_sum = utility_sum + piecewise_utility(x_tilde(i,t), ...
                D(i,t), beta, sigma);
        end

        % 系统成本
        grid_cost_t = C_grid_func(G(t), t);
        battery_cost_t = C_battery_func(B(t), SOC(t), battery_params);

        total_cost = total_cost - utility_sum + grid_cost_t + battery_cost_t;

        cost_history(k) = total_cost;

        %% 收敛检查
        eps_primal = sqrt(T) * abs_tol + rel_tol * ...
            max(norm(sum(x_tilde, 1), 2), norm(z, 2));
        eps_dual = sqrt(T) * abs_tol + rel_tol * norm(lambda, 2);

        if primal_residual(k) < eps_primal && dual_residual(k) < eps_dual
            fprintf('在第 %d 次迭代收敛\n', k);
            break;
        end

        z_prev = z;

        %% 自适应调整rho
        if mod(k, 50) == 0
            if primal_residual(k) > 10 * dual_residual(k)
                rho = rho * 1.5;
                fprintf('迭代 %d: 增加rho到 %.2f\n', k, rho);
            elseif dual_residual(k) > 10 * primal_residual(k)
                rho = rho / 1.5;
                fprintf('迭代 %d: 减少rho到 %.2f\n', k, rho);
            end
        end

        if mod(k, 100) == 0
            fprintf('迭代 %d: 成本=%.4f, 原始残差=%.4e\n', ...
                k, cost_history(k), primal_residual(k));
        end
    end
end

%% 最终结果
x = x_tilde;
fprintf('\n=== 优化结果汇总 ===\n');
fprintf('总成本: %.4f ¥\n', cost_history(k));
fprintf('总电网用电: %.2f kWh\n', sum(G));
fprintf('总电池循环: %.2f kWh\n', sum(abs(B)));
fprintf('平均SOC: %.2f\n', mean(SOC));

%% 可视化
figure('Position', [100, 100, 1400, 800]);

% 1. 成本分解
subplot(2, 3, 1);
cost_utility = 0;
cost_grid = 0;
cost_battery = 0;
for t = 1:T
    for i = 1:n
        cost_utility = cost_utility - piecewise_utility(x(i,t), D(i,t), beta, sigma);
    end
    cost_grid = cost_grid + C_grid_func(G(t), t);
    cost_battery = cost_battery + C_battery_func(B(t), SOC(t), battery_params);
end
pie([-cost_utility, cost_grid, cost_battery]);
title('成本分解');
legend({'用户效用成本', '电网成本', '电池成本'}, 'Location', 'best');

% 2. 电力平衡
subplot(2, 3, 2);
plot(1:T, sum(x,1), 'b-', 'LineWidth', 2); hold on;
plot(1:T, PV+G+max(0,B)-min(0,B), 'r--', 'LineWidth', 2);
xlabel('时间 (小时)');
ylabel('功率 (kW)');
title('电力平衡');
legend({'总需求', '总供应'}, 'Location', 'best');
grid on;

% 3. 电源组成
subplot(2, 3, 3);
area(1:T, [PV; G; max(0,B); min(0,B)]');
xlabel('时间 (小时)');
ylabel('功率 (kW)');
title('电源组成');
legend({'光伏', '电网', '电池放电', '电池充电'});
grid on;

% 4. SOC和电价
subplot(2, 3, 4);
yyaxis left;
plot(1:T, SOC, 'b-', 'LineWidth', 2);
ylabel('SOC');
ylim([0 1]);
yyaxis right;
plot(1:T, time_of_use_price, 'r-', 'LineWidth', 2);
ylabel('电价 (¥/kWh)');
xlabel('时间 (小时)');
title('SOC和分时电价');
grid on;

% 5. 收敛过程
subplot(2, 3, 5);
semilogy(1:k, primal_residual(1:k), 'b-*', 'LineWidth', 2); hold on;
semilogy(1:k, dual_residual(1:k), 'r--*', 'LineWidth', 2);
xlabel('迭代次数');
ylabel('残差');
title('ADMM收敛过程');
legend({'原始残差', '对偶残差'});

% 6. 用户用电对比
subplot(2, 3, 6);
user_to_plot = min(5, n);
for i = 1:user_to_plot
    plot(1:T, D(i,1:T), '-', 'LineWidth', 1); hold on;
    plot(1:T, x(i,1:T), '--', 'LineWidth', 1.5);
end
xlabel('时间 (小时)');
ylabel('用电量 (kW)');
title(sprintf('前%d个用户用电对比', user_to_plot));
grid on;

%% 辅助函数
function U = piecewise_utility(x, D, beta, sigma)
if x <= D
    U = beta * log(1 + sigma * x);
else
    U = beta * log(1 + sigma * D);
end
end

function cost = user_subproblem(xi, Di, lambda_i, z_i, beta, sigma, rho, n)
U = piecewise_utility(xi, Di, beta, sigma);
consensus_term = lambda_i * (xi/n) + (rho/2) * (xi/n - z_i/n)^2;
cost = -U + consensus_term;
end

function cost = battery_total_cost(B, SOC, params)
% 电池总成本函数
cycle_cost = params.cycle_cost_coeff * abs(B);
SOC_penalty = params.SOC_penalty * (SOC - params.SOC_target)^2;

if B > 0 % 放电
    efficiency_cost = params.discharge_loss * B;
else % 充电
    efficiency_cost = params.charge_loss * abs(B);
end

% 功率超限惩罚
power_penalty = 0;
if abs(B) > params.B_max
    power_penalty = 1000 * (abs(B) - params.B_max)^2;
end

cost = cycle_cost + SOC_penalty + efficiency_cost + power_penalty;
end

function cost = system_cost(G, B, SOC, t, C_grid, C_battery, params)
% 系统总成本
grid_cost = C_grid(G, t);
battery_cost = C_battery(B, SOC, params);
cost = grid_cost + battery_cost;
end

function [c, ceq] = power_balance_constraint(vars, net_demand)
% 功率平衡约束
G = vars(1);
B = vars(2);
ceq = G + B - net_demand;  % G + B = net_demand
c = [];
end