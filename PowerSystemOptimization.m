classdef PowerSystemOptimization
    % PowerSystemOptimization 类：电力系统优化调度主程序
    
    properties
        % 基本参数
        n           % 用户数量
        T           % 时间周期数（24小时）
        
        % 用户效用参数
        beta
        sigma
        
        % 电网成本参数
        a0
        b0
        c0
        
        % 电池参数
        E_max
        B_max
        SOC_min
        SOC_max
        SOC_initial
        
        % ADMM参数
        rho
        max_admm_iter
        abs_tol
        rel_tol
        
        % 模拟数据
        D           % 用户基准需求
        PV          % 光伏发电预测
        time_price  % 分时电价
        
        % 优化结果
        X_opt       % 最优用电量 [kW]
        G_opt       % 最优电网供电 [kW]
        B_opt       % 最优电池充放电 [kW]
        SOC_opt     % SOC轨迹
        Lambda_opt  % 最优对偶变量
        convergence_info % 收敛信息
        total_cost_all   % 总成本
        
        % 每个时段的优化问题对象
        hourly_optimizers
    end
    
    methods
        function obj = PowerSystemOptimization()
            % 构造函数：初始化参数和生成模拟数据
            
            %% 参数设置
            obj.n = 30;           % 用户数量
            obj.T = 24;          % 时间周期数（24小时）
            obj.beta = 1.5;      % 效用函数系数
            obj.sigma = 0.8;     % 效用函数形状参数
            
            % 电网成本参数
            obj.a0 = 0.001;       % 二次项系数 (¥/kW²)
            obj.b0 = 0.02;        % 一次项系数 (¥/kW)  
            obj.c0 = 0;          % 常数项
            
            % 电池参数
            obj.E_max = 1000;     % 电池容量 (kWh)
            obj.B_max = 250;      % 最大充放电功率 (kW)
            obj.SOC_min = 0.15;  % 最小SOC
            obj.SOC_max = 0.95;  % 最大SOC
            obj.SOC_initial = 0.5; % 初始SOC
            
            % ADMM参数
            obj.rho = 0.8;       % 初始惩罚参数
            obj.max_admm_iter = 200; % 每个时段最大ADMM迭代次数
            obj.abs_tol = 1e-5;
            obj.rel_tol = 1e-5;
            
            %% 生成模拟数据
            obj = obj.generate_simulated_data();
            
            %% 初始化存储结果
            obj.X_opt = zeros(obj.n, obj.T);
            obj.G_opt = zeros(1, obj.T);
            obj.B_opt = zeros(1, obj.T);
            obj.SOC_opt = zeros(1, obj.T);
            obj.Lambda_opt = zeros(1, obj.T);
            obj.convergence_info = zeros(obj.T, 2);
            obj.total_cost_all = 0;
            obj.hourly_optimizers = cell(1, obj.T);
        end
        
        function obj = generate_simulated_data(obj)
            % 生成模拟数据
            
            rng(2023); % 设置随机种子
            
            % 生成用户基准需求 D^A_{i,t}
            obj.D = zeros(obj.n, obj.T);
            for i = 1:obj.n
                % 每个用户有不同的用电模式
                user_type = randi([1, 3]); % 1:住宅, 2:商业, 3:工业
                base_load = [1.5, 3.0, 8.0]; % 基础负荷
                
                % 生成日负荷曲线
                t_vec = 1:obj.T;
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
                
                obj.D(i, :) = base_load(user_type) * pattern + 0.2*randn(1, obj.T);
                obj.D(i, :) = max(0.1, obj.D(i, :)); % 确保非负
            end
            
            % 光伏发电预测 PV_A
            obj.PV = zeros(1, obj.T);
            daylight_hours = 6:19; % 白天时段
            sun_power = 25*7; % 峰值光伏功率
            obj.PV(daylight_hours) = sun_power * sin(pi*(daylight_hours-6)/(19-6)).^2;
            obj.PV = obj.PV + 2*randn(1, obj.T); % 添加噪声
            obj.PV = max(0, obj.PV); % 确保非负
            
            % 分时电价
            obj.time_price = zeros(1, obj.T);
            peak_hours = [8:11, 18:21];    % 峰时
            flat_hours = [7, 12:17, 22];   % 平时
            valley_hours = [1:6, 23:24];   % 谷时
            obj.time_price(peak_hours) = 0.8;  % 峰时电价
            obj.time_price(flat_hours) = 0.5;  % 平时电价
            obj.time_price(valley_hours) = 0.3; % 谷时电价
        end
        
        function obj = run_optimization(obj)
            % 运行优化调度
            
            %% 主循环：按时间顺序处理
            SOC_current = obj.SOC_initial;
            
            for t = 1:obj.T
                % 创建当前时段的优化问题对象
                optimizer = HourlyOptimization(...
                    t, obj.D(:, t), obj.PV(t), obj.time_price(t), SOC_current, ...
                    obj.beta, obj.sigma, obj.a0, obj.b0, obj.c0, ...
                    obj.B_max, obj.E_max, obj.SOC_min, obj.SOC_max, ...
                    obj.rho, obj.max_admm_iter, obj.abs_tol, obj.rel_tol);
                
                % 求解当前时段的优化问题
                optimizer = optimizer.solve();
                
                % 存储优化结果
                obj.X_opt(:, t) = optimizer.x;
                obj.G_opt(t) = optimizer.G;
                obj.B_opt(t) = optimizer.B;
                obj.SOC_opt(t) = optimizer.SOC_new;
                obj.Lambda_opt(t) = optimizer.lambda;
                obj.convergence_info(t, 1) = optimizer.iter_used;
                obj.convergence_info(t, 2) = optimizer.primal_res_final;
                
                % 存储优化器对象
                obj.hourly_optimizers{t} = optimizer;
                
                % 更新SOC供下一个时段使用
                SOC_current = optimizer.SOC_new;
                
                % 计算当前时段成本并累加
                [utility_val, grid_cost_val, battery_cost_val] = ...
                    obj.calculate_hourly_costs(t);
                total_cost_t = -utility_val + grid_cost_val + battery_cost_val;
                obj.total_cost_all = obj.total_cost_all + total_cost_t;
            end
            
            %% 计算总体性能指标
            obj.display_results();
        end
        
        function [utility_val, grid_cost_val, battery_cost_val] = ...
                calculate_hourly_costs(obj, t)
            % 计算当前时段的成本
            
            % 用户效用
            utility_val = 0;
            for i = 1:obj.n
                utility_val = utility_val + obj.piecewise_log_utility(...
                    obj.X_opt(i, t), obj.D(i, t));
            end
            
            % 电网成本
            grid_cost_val = obj.a0*obj.G_opt(t)^2 + ...
                (obj.b0 + obj.time_price(t))*obj.G_opt(t) + obj.c0;
            
            % 电池成本
            battery_cost_val = obj.battery_cost(obj.B_opt(t), obj.SOC_opt(t));
        end
        
        function U = piecewise_log_utility(obj, x, D)
            % 分段对数效用函数
            if x <= D
                U = obj.beta * log(1 + obj.sigma * x);
            else
                U = obj.beta * log(1 + obj.sigma * D);
            end
        end
        
        function cost = battery_cost(obj, x, SOC)
            % 电池成本函数
            if x <= 0  % 充电
                cost = (1-SOC)/SOC * exp(0) * 0.1 * (exp(x) - 1);
            else  % 放电
                cost = 0.00011*x^2 + (1-SOC)/SOC * exp(0) * 0.1 * x;
            end
        end
        
        function display_results(obj)
            % 显示优化结果
            
            fprintf('\n\n=== 优化完成 ===\n');
            fprintf('总成本: %.2f ¥\n', obj.total_cost_all);
            fprintf('平均每时段成本: %.2f ¥\n', obj.total_cost_all/obj.T);
            fprintf('总用电量: %.2f kWh\n', sum(sum(obj.X_opt)));
            fprintf('总电网供电: %.2f kWh\n', sum(obj.G_opt));
            fprintf('总电池充电: %.2f kWh\n', sum(obj.B_opt(obj.B_opt < 0)));
            fprintf('总电池放电: %.2f kWh\n', sum(obj.B_opt(obj.B_opt > 0)));
            fprintf('最终SOC: %.3f\n', obj.SOC_opt(end));
            fprintf('平均ADMM迭代次数: %.1f\n', mean(obj.convergence_info(:,1)));
            
            %% 成本分解分析
            [total_utility, total_grid_cost, total_battery_cost] = ...
                obj.analyze_costs();
            
            fprintf('\n=== 成本分解 ===\n');
            fprintf('用户效用总价值: %.2f ¥\n', total_utility);
            fprintf('电网总成本: %.2f ¥ (%.1f%%)\n', ...
                total_grid_cost, 100*total_grid_cost/obj.total_cost_all);
            fprintf('电池总成本: %.2f ¥ (%.1f%%)\n', ...
                total_battery_cost, 100*total_battery_cost/obj.total_cost_all);
            fprintf('净总成本: %.2f ¥\n', ...
                -total_utility + total_grid_cost + total_battery_cost);
        end
        
        function [total_utility, total_grid_cost, total_battery_cost] = ...
                analyze_costs(obj)
            % 分析总成本分解
            
            total_utility = 0;
            total_grid_cost = 0;
            total_battery_cost = 0;
            
            for t = 1:obj.T
                % 用户效用
                for i = 1:obj.n
                    total_utility = total_utility + ...
                        obj.piecewise_log_utility(obj.X_opt(i,t), obj.D(i,t));
                end
                
                % 电网成本
                total_grid_cost = total_grid_cost + ...
                    obj.a0*obj.G_opt(t)^2 + (obj.b0 + obj.time_price(t))*obj.G_opt(t) + obj.c0;
                
                % 电池成本
                total_battery_cost = total_battery_cost + ...
                    obj.battery_cost(obj.B_opt(t), obj.SOC_opt(t));
            end
        end
        
        function visualize_results(obj)
            % 可视化所有结果
            
            visualize_results_impl(obj.X_opt, obj.G_opt, obj.B_opt, ...
                obj.PV, obj.SOC_opt, obj.D, obj.time_price, ...
                obj.convergence_info, obj.total_cost_all, obj.Lambda_opt);
        end
        
        function save_results(obj, filename)
            % 保存结果到文件
            
            if nargin < 2
                filename = 'optimization_results.mat';
            end
            
            save(filename, 'obj');
            fprintf('结果已保存到 %s\n', filename);
        end
    end
end