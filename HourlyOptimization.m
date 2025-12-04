classdef HourlyOptimization
    % HourlyOptimization 类：封装每个小时的优化问题
    
    properties
        % 时段索引
        t
        % 当前时段数据
        D_t
        PV_t
        time_price_t
        SOC_current
        
        % 优化参数
        beta
        sigma
        a0
        b0
        c0
        B_max
        E_max
        SOC_min
        SOC_max
        
        % ADMM参数
        rho
        max_admm_iter
        abs_tol
        rel_tol
        
        % 优化结果
        x
        G
        B
        lambda
        SOC_new
        iter_used
        primal_res_final
        converged
        cost_history
        primal_res_history
        dual_res_history
    end
    
    methods
        function obj = HourlyOptimization(t, D_t, PV_t, time_price_t, SOC_current, ...
                beta, sigma, a0, b0, c0, B_max, E_max, SOC_min, SOC_max, ...
                rho, max_admm_iter, abs_tol, rel_tol)
            % 构造函数：初始化时段优化问题
            
            obj.t = t;
            obj.D_t = D_t;
            obj.PV_t = PV_t;
            obj.time_price_t = time_price_t;
            obj.SOC_current = SOC_current;
            
            obj.beta = beta;
            obj.sigma = sigma;
            obj.a0 = a0;
            obj.b0 = b0;
            obj.c0 = c0;
            obj.B_max = B_max;
            obj.E_max = E_max;
            obj.SOC_min = SOC_min;
            obj.SOC_max = SOC_max;
            
            obj.rho = rho;
            obj.max_admm_iter = max_admm_iter;
            obj.abs_tol = abs_tol;
            obj.rel_tol = rel_tol;
            
            % 初始化优化结果
            obj.x = D_t;
            obj.G = max(0, sum(obj.x) - PV_t);
            obj.B = 0;
            obj.lambda = 0;
            obj.SOC_new = SOC_current;
            obj.iter_used = max_admm_iter;
            obj.primal_res_final = 0;
            obj.converged = false;
        end
        
        function obj = solve(obj)
            % 求解当前时段的优化问题
            
            fprintf('\n=== 时段 %02d:00 - %02d:00 ===\n', obj.t-1, obj.t);
            fprintf('光伏预测: %.2f kW, 分时电价: %.2f ¥/kWh\n', ...
                obj.PV_t, obj.time_price_t);
            
            %% ADMM初始化（当前时段）
            rho_t = obj.rho;  % 当前时段的惩罚参数
            x = obj.D_t;      % 初始化为基准需求
            G = max(0, sum(x) - obj.PV_t); % 初始电网供电
            B = 0;            % 初始电池功率 
            lambda = 0;       % 对偶变量
            z_prev = G + obj.PV_t + B;
            
            %% 存储迭代历史
            primal_res_history = zeros(obj.max_admm_iter, 1);
            dual_res_history = zeros(obj.max_admm_iter, 1);
            cost_history = zeros(obj.max_admm_iter, 1);
            
            %% ADMM迭代（当前时段）
            converged = false;
            iter_used = obj.max_admm_iter;
            
            for k = 1:obj.max_admm_iter
                % 保存上一次的z值（用于计算对偶残差）
                if k > 1
                    z_prev = z;
                end
                
                %% 步骤1: 更新x（用户侧）
                total_x_prev = sum(x);
                
                % 并行更新每个用户的用电量
                n = length(obj.D_t);
                for i = 1:n  
                    % 用户i的本地优化问题
                    z0 = sum(x) - x(i) - (G + obj.PV_t + B);
                    fun = @(xi) obj.user_local_optimization(...
                        xi, obj.D_t(i), lambda, z0);
                    
                    % 边界约束
                    lb = 0;
                    ub = 1 * obj.D_t(i); % 最大不超过2.5倍基准需求
                    
                    % 求解一维优化问题
                    options = optimset('Display', 'off', 'TolX', 1e-6);
                    x(i) = fminbnd(fun, lb, ub, options);
                end
                
                %% 步骤2: 更新系统变量 (z, G, B) 
                
                % 系统优化：最小化电网成本+电池成本+ADMM惩罚项
                [G_new, B_new, feasible, ~] = obj.system_optimization(...
                    lambda, sum(x), rho_t);
                
                if feasible
                    G = G_new;
                    B = B_new; 
                    z = G + obj.PV_t + B;  
                else
                    warning('时段 %d 迭代 %d: 系统问题不可行，保持上一解', obj.t, k);
                end
                
                %% 步骤3: 更新对偶变量
                primal_res = sum(x) - (G + obj.PV_t + B);
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
                    user_utility = user_utility + obj.piecewise_log_utility(...
                        x(i), obj.D_t(i));
                end
                
                grid_cost_val = obj.C_grid_func(G);
                battery_cost_val = obj.battery_cost(B, obj.SOC_current); 
                
                cost_history(k) = -user_utility + grid_cost_val + ...
                    battery_cost_val;
                
                %% 收敛检查
                eps_primal = obj.abs_tol + obj.rel_tol * max(abs(sum(x)), ...
                    abs(G + obj.PV_t + B));
                eps_dual = obj.abs_tol + obj.rel_tol * abs(z);
                
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
                if mod(k, 20) == 0 && k < obj.max_admm_iter
                    fprintf('  迭代 %3d: 原始残差=%.2e, 对偶残差=%.2e, 成本=%.2f\n', ...
                        k, primal_res_history(k), dual_res_history(k), cost_history(k));
                end
            end
            
            if ~converged
                fprintf('  ADMM未收敛，使用最后迭代结果\n');
            end
            
            %% 存储优化结果
            obj.x = x;
            obj.G = G;
            obj.B = B;
            obj.lambda = lambda;
            obj.iter_used = iter_used;
            obj.primal_res_final = primal_res_history(min(iter_used, obj.max_admm_iter));
            obj.converged = converged;
            obj.cost_history = cost_history(1:iter_used);
            obj.primal_res_history = primal_res_history(1:iter_used);
            obj.dual_res_history = dual_res_history(1:iter_used);
            
            % 更新SOC（供下一个时段使用）
            obj.SOC_new = obj.update_battery_SOC(obj.SOC_current, B);
            
            % 计算当前时段成本
            user_utility_t = 0;
            n = length(obj.D_t);
            for i = 1:n
                user_utility_t = user_utility_t + obj.piecewise_log_utility(...
                    x(i), obj.D_t(i));
            end
            
            grid_cost_t = obj.C_grid_func(G);
            battery_cost_t = obj.battery_cost(B, obj.SOC_current);
            total_cost_t = -user_utility_t + grid_cost_t + battery_cost_t;
            
            fprintf('  优化结果: 总需求=%.2f kW, 电网=%.2f kW, 电池=%.2f kW, SOC=%.3f\n', ...
                sum(x), G, B, obj.SOC_new);
            fprintf('  时段成本: 用户效用=%.2f, 电网成本=%.2f, 电池成本=%.2f, 总成本=%.2f\n', ...
                user_utility_t, grid_cost_t, battery_cost_t, total_cost_t);
        end
        
        function cost = battery_cost(obj, x, SOC)
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
        
        function U = piecewise_log_utility(obj, x, D)
            % 分段对数效用函数
            if x <= D
                U = obj.beta * log(1 + obj.sigma * x);
            else
                U = obj.beta * log(1 + obj.sigma * D);
            end
        end
        
        function cost = user_local_optimization(obj, xi, Di, lambda, z)
            % 用户本地优化目标函数
            utility = obj.piecewise_log_utility(xi, Di);
            consensus_term = lambda * (xi) + (obj.rho/2) * (xi + z)^2;
            cost = -utility + consensus_term;
        end
        
        function cost = C_grid_func(obj, x)
            % 电网成本函数（考虑分时电价）
            cost = obj.a0*x^2 + (obj.b0 + obj.time_price_t)*x + obj.c0;
        end
        
        function [G_opt, B_opt, feasible, total_cost] = system_optimization(...
                obj, lambda, sum_x, rho_t)
            % 系统优化：求解最优的G和B
            
            % 净需求
            net_demand = sum_x - obj.PV_t;
            
            % 定义优化问题
            options = optimoptions('fmincon', 'Display', 'off', ...
                'Algorithm', 'sqp', 'MaxIterations', 100);
            
            % 初始猜测
            x0 = [max(0, net_demand), 0];  % [G, B]
            
            % 边界约束
            lb = [0, -obj.B_max];
            ub = [Inf, obj.B_max];
            
            % 线性等式约束：G + B = net_demand
            Aeq = [1, 1];
            beq = net_demand;
            
            % 目标函数
            function f = system_objective(vars)
                G = vars(1);
                B = vars(2);
                
                % 电网成本（考虑分时电价）
                grid_cost = obj.a0*G^2 + (obj.b0 + obj.time_price_t)*G;
                
                % 电池成本
                battery_cost_val = obj.battery_cost(B, obj.SOC_current);
                
                % ADMM惩罚项  
                admm_penalty = lambda * (-G - obj.PV_t - B) + ...
                    (rho_t/2) * (sum_x - (G + obj.PV_t + B))^2;
                
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
        
        function SOC_new = update_battery_SOC(obj, SOC_old, B)
            % 更新电池SOC
            % B: 电池功率 (kW), 正为放电，负为充电
            % 假设时间间隔为1小时
            
            eta_ch = 0.95;  % 充电效率
            eta_dis = 0.97; % 放电效率
            
            if B > 0  % 放电
                delta_SOC = -B / (eta_dis * obj.E_max);
            else  % 充电
                delta_SOC = -B * eta_ch / obj.E_max;
            end
            
            SOC_new = SOC_old + delta_SOC;
            SOC_new = max(obj.SOC_min, min(obj.SOC_max, SOC_new));
        end
    end
end