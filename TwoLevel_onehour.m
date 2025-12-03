classdef TwoLevel_onehour < handle

    properties
        % 系统模型
        MA; MB; MC;

        % 优化变量
        tDA; tDB; tDC;       % 需求变量
        tSA; tSB; tSC;       % 供给变量 [G, T1, T2, Eb]

        % ========== 外层ADMM变量（功率平衡） ==========
        % 对偶变量（价格）
        lambda;              % [λ1, λ2, λ3] - 功率平衡对偶变量
        
        % 共识变量 - 修正：每个区域应该有单独的需求和供给共识变量
        d_consensus;         % 需求共识变量 [D1_cons, D2_cons, D3_cons]
        s_consensus;         % 供给共识变量 [S1_cons, S2_cons, S3_cons]
        
        % 对偶变量（用于共识约束）
        u_demand;            % 需求共识的对偶变量
        u_supply;            % 供给共识的对偶变量

        % 残差和参数
        rho_outer = 2.0;     % 外层惩罚参数 - 增加初始值
        eps_outer = 1e-3;    % 外层收敛容差 - 放宽要求
        max_iter_outer = 100;% 增加外层迭代次数

        % ========== 内层ADMM变量（传输一致性） ==========
        % 传输功率的全局共识变量
        z_ab; z_ac; z_bc;    % 全局传输功率

        % 对偶变量
        mu_ab; mu_ac; mu_bc; % 传输一致性对偶变量

        % 残差和参数
        rho_inner = 2.0;     % 内层惩罚参数 - 增加初始值
        eps_inner = 1e-3;    % 内层收敛容差 - 放宽要求
        max_iter_inner = 100;% 增加内层迭代次数

        % 结果记录
        TotalCost;
        outer_iter_count = 0;
        total_inner_iter = 0;

        % 成本参数
        CofT = 0.5;
        
        % 收敛历史记录（用于调试）
        outer_residuals = [];
        inner_residuals = [];
        prices_history = [];
    end

    methods
        function obj = TwoLevel_onehour(NT)
            % 构造函数
            [MA, MB, MC] = date_generation(NT);
            obj.MA = MA;
            obj.MB = MB;
            obj.MC = MC;

            % 初始化变量
            obj.lambda = ones(1, 3) * 0.5;  % 初始价格设为0.5
            
            % 初始化共识变量
            obj.d_consensus = zeros(1, 3);
            obj.s_consensus = zeros(1, 3);
            obj.u_demand = zeros(1, 3);
            obj.u_supply = zeros(1, 3);

            obj.tDA = zeros(obj.MA.N, 1);
            obj.tDB = zeros(obj.MB.N, 1);
            obj.tDC = 0;
            obj.tSA = [0.1, 0, 0, 0];  % 更好的初始值
            obj.tSB = [0.1, 0, 0, 0];
            obj.tSC = [0.1, 0, 0, 0];

            % 初始化传输变量
            obj.z_ab = 0; obj.z_ac = 0; obj.z_bc = 0;
            obj.mu_ab = 0; obj.mu_ac = 0; obj.mu_bc = 0;
        end

        function Solve(obj, t)
            % 双层ADMM主函数
            obj.outer_iter_count = 0;
            obj.total_inner_iter = 0;
            obj.outer_residuals = [];
            obj.inner_residuals = [];
            obj.prices_history = [];

            for outer_iter = 1:obj.max_iter_outer
                obj.outer_iter_count = outer_iter;

                fprintf('\n=== Outer ADMM Iteration %d (ρ_outer=%.3f, ρ_inner=%.3f) ===\n', ...
                    outer_iter, obj.rho_outer, obj.rho_inner);

                % ========== Step 1: 更新需求（修正的ADMM形式） ==========
                obj.UpdateDemand_Corrected(t);

                % ========== Step 2: 更新供给（内层ADMM） ==========
                inner_iters = obj.InnerADMM_Corrected(t);
                obj.total_inner_iter = obj.total_inner_iter + inner_iters;

                % ========== Step 3: 更新外层共识变量 ==========
                obj.UpdateOuterConsensus_Corrected(t);

                % ========== Step 4: 更新外层对偶变量（价格和共识对偶变量） ==========
                obj.UpdateOuterDualVariables_Corrected(t);

                % ========== Step 5: 计算并记录残差 ==========
                [outer_r_norm, outer_s_norm, power_balance] = obj.ComputeOuterResiduals(t);
                obj.outer_residuals = [obj.outer_residuals; outer_r_norm, outer_s_norm];
                obj.prices_history = [obj.prices_history; obj.lambda];

                % 显示当前状态
                fprintf('Outer residuals: r=%.4e, s=%.4e, power imbalance=%.4e\n', ...
                    outer_r_norm, outer_s_norm, norm(power_balance));
                fprintf('Prices: λ=[%.4f, %.4f, %.4f]\n', obj.lambda(1), obj.lambda(2), obj.lambda(3));

                % ========== Step 6: 检查外层收敛 ==========
                if obj.CheckOuterConvergence_Corrected(outer_r_norm, outer_s_norm, power_balance)
                    fprintf('✓ Outer ADMM converged in %d iterations\n', outer_iter);
                    break;
                end

                % ========== Step 7: 自适应调整惩罚参数 ==========
                obj.AdaptParameters_Corrected(outer_iter, outer_r_norm, outer_s_norm);
                
                % 可选：绘图查看当前状态
                if mod(outer_iter, 10) == 0
                    obj.draw(t);
                    pause(0.1);
                end
            end

            % 计算总成本
            obj.CalculateTotalCost(t);

            % 显示结果
            obj.DisplayResults_Corrected(t);
            
            % 绘制收敛曲线
            obj.PlotConvergence();
        end

        function UpdateDemand_Corrected(obj, t)
            % 修正的需求更新函数
            % ADMM形式：min -U(x) + λ*x + (ρ/2)||x - d_consensus + u_demand||²
            
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

            % 区域A的需求
            for i = 1:obj.MA.N
                b = obj.MA.D_it(i, t);
                
                % 修正的目标函数
                fun = @(x) -obj.U_fun(x, b) + obj.lambda(1)*x ...
                          + (obj.rho_outer/2) * (x - obj.d_consensus(1) + obj.u_demand(1))^2;
                
                % 使用更好的初始值
                x0 = min(max(obj.tDA(i), 0), b);
                [obj.tDA(i), ~] = fmincon(fun, x0, [], [], [], [], 0, b, [], options);
            end

            % 区域B的需求
            for j = 1:obj.MB.N
                b = obj.MB.D_it(j, t);
                
                fun = @(x) -obj.U_fun(x, b) + obj.lambda(2)*x ...
                          + (obj.rho_outer/2) * (x - obj.d_consensus(2) + obj.u_demand(2))^2;
                
                x0 = min(max(obj.tDB(j), 0), b);
                [obj.tDB(j), ~] = fmincon(fun, x0, [], [], [], [], 0, b, [], options);
            end

            % 区域C的需求
            b = obj.MC.D_it(1, t);
            fun = @(x) -obj.U_fun(x, b) + obj.lambda(3)*x ...
                      + (obj.rho_outer/2) * (x - obj.d_consensus(3) + obj.u_demand(3))^2;
            
            x0 = min(max(obj.tDC, 0), b);
            [obj.tDC, ~] = fmincon(fun, x0, [], [], [], [], 0, b, [], options);
        end

        function inner_iterations = InnerADMM_Corrected(obj, t)
            % 修正的内层ADMM
            
            inner_iterations = 0;
            z_old = [obj.z_ab; obj.z_ac; obj.z_bc];
            
            for inner_iter = 1:obj.max_iter_inner
                inner_iterations = inner_iter;
                z_old_iter = [obj.z_ab; obj.z_ac; obj.z_bc];

                % ===== Step 1: 更新各区域供给变量 =====
                obj.UpdateRegionsSupply_Corrected(t);

                % ===== Step 2: 更新全局共识变量 =====
                obj.UpdateGlobalTransmission_Corrected();

                % ===== Step 3: 更新内层对偶变量 =====
                obj.UpdateInnerDualVariables_Corrected();

                % ===== Step 4: 计算残差 =====
                [r_norm, s_norm] = obj.ComputeInnerResiduals_Corrected(z_old_iter);
                
                if inner_iter == 1
                    obj.inner_residuals = [obj.inner_residuals; r_norm, s_norm];
                end

                % ===== Step 5: 检查收敛 =====
                if r_norm < obj.eps_inner && s_norm < obj.eps_inner
                    if inner_iter > 1
                        fprintf('  Inner ADMM converged in %d iterations (r=%.2e, s=%.2e)\n', ...
                            inner_iter, r_norm, s_norm);
                    end
                    break;
                end
                
                % ===== Step 6: 自适应调整内层参数 =====
                if mod(inner_iter, 5) == 0
                    obj.AdaptInnerRho_Corrected(r_norm, s_norm);
                end
            end
        end

        function UpdateRegionsSupply_Corrected(obj, t)
            % 修正的供给更新函数
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

            % 区域A - 修正目标函数
            funA = @(x) obj.RegionA_Objective_Corrected(x, t);
            lbA = [0, -obj.MA.Tmax, -obj.MA.Tmax, -obj.MA.BC*(1-obj.MA.SOC(t))];
            ubA = [obj.MA.Gmax, obj.MA.Tmax, obj.MA.Tmax, obj.MA.BC*obj.MA.SOC(t)];
            [obj.tSA, ~] = fmincon(funA, obj.tSA, [], [], [], [], lbA, ubA, [], options);

            % 区域B
            funB = @(x) obj.RegionB_Objective_Corrected(x, t);
            lbB = [0, -obj.MB.Tmax, -obj.MB.Tmax, -obj.MB.BC*(1-obj.MB.SOC(t))];
            ubB = [obj.MB.Gmax, obj.MB.Tmax, obj.MB.Tmax, obj.MB.BC*obj.MB.SOC(t)];
            [obj.tSB, ~] = fmincon(funB, obj.tSB, [], [], [], [], lbB, ubB, [], options);

            % 区域C
            funC = @(x) obj.RegionC_Objective_Corrected(x, t);
            lbC = [0, -obj.MC.Tmax, -obj.MC.Tmax, -obj.MC.BC*(1-obj.MC.SOC(t))];
            ubC = [obj.MC.Gmax, obj.MC.Tmax, obj.MC.Tmax, obj.MC.BC*obj.MC.SOC(t)];
            [obj.tSC, ~] = fmincon(funC, obj.tSC, [], [], [], [], lbC, ubC, [], options);
        end

        function f = RegionA_Objective_Corrected(obj, x, t)
            % 修正的区域A目标函数
            G = x(1); T_ab = x(2); T_ac = x(3); Eb = x(4);
            
            % 计算实际供给
            S_actual = G + obj.MA.PV(t) - T_ab - T_ac + Eb;
            
            % 基本成本项
            f = obj.G_fun(G) ...
                + obj.E_fun(Eb, obj.MA.SOC(t), -obj.MA.BC*(1-obj.MA.SOC(t))) ...
                + obj.CofT*(T_ab^2 + T_ac^2);
            
            % 价格项（注意符号）
            f = f - obj.lambda(1) * S_actual;
            
            % 外层ADMM惩罚项（供给共识）
            f = f + (obj.rho_outer/2) * (S_actual - obj.s_consensus(1) + obj.u_supply(1))^2;
            
            % 内层ADMM惩罚项（传输一致性）
            f = f + (obj.rho_inner/2) * ((T_ab - obj.z_ab + obj.mu_ab)^2 ...
                                       + (T_ac - obj.z_ac + obj.mu_ac)^2);
        end

        function f = RegionB_Objective_Corrected(obj, x, t)
            % 修正的区域B目标函数
            G = x(1); T_ab_recv = x(2); T_bc = x(3); Eb = x(4);
            
            % 注意：T_ab_recv是接收的功率（应为负值或正值？需要统一符号）
            % 假设T_ab_recv为接收功率的绝对值，所以是正值
            S_actual = G + obj.MB.PV(t) + T_ab_recv - T_bc + Eb;
            
            f = obj.G_fun(G) ...
                + obj.E_fun(Eb, obj.MB.SOC(t), -obj.MB.BC*(1-obj.MB.SOC(t))) ...
                + obj.CofT*(T_ab_recv^2 + T_bc^2);
            
            f = f - obj.lambda(2) * S_actual;
            
            % 外层ADMM惩罚项
            f = f + (obj.rho_outer/2) * (S_actual - obj.s_consensus(2) + obj.u_supply(2))^2;
            
            % 内层ADMM惩罚项
            % 注意：接收功率应该取负号，因为传输方向相反
            f = f + (obj.rho_inner/2) * ((-T_ab_recv - obj.z_ab + obj.mu_ab)^2 ...
                                       + (T_bc - obj.z_bc + obj.mu_bc)^2);
        end

        function f = RegionC_Objective_Corrected(obj, x, t)
            % 修正的区域C目标函数
            G = x(1); T_ac_recv = x(2); T_bc_recv = x(3); Eb = x(4);
            
            S_actual = G + obj.MC.PV(t) + T_ac_recv + T_bc_recv + Eb;
            
            f = obj.G_fun(G) ...
                + obj.E_fun(Eb, obj.MC.SOC(t), -obj.MC.BC*(1-obj.MC.SOC(t))) ...
                + obj.CofT*(T_ac_recv^2 + T_bc_recv^2);
            
            f = f - obj.lambda(3) * S_actual;
            
            % 外层ADMM惩罚项
            f = f + (obj.rho_outer/2) * (S_actual - obj.s_consensus(3) + obj.u_supply(3))^2;
            
            % 内层ADMM惩罚项
            f = f + (obj.rho_inner/2) * ((-T_ac_recv - obj.z_ac + obj.mu_ac)^2 ...
                                       + (-T_bc_recv - obj.z_bc + obj.mu_bc)^2);
        end

        function UpdateGlobalTransmission_Corrected(obj)
            % 修正的全局传输共识变量更新
            % 使用加权平均而不是简单平均
            
            alpha = 0.5;  % 加权系数
            
            % 对于A->B
            obj.z_ab = alpha * obj.z_ab + (1-alpha) * 0.5 * (obj.tSA(2) - obj.tSB(2) + obj.mu_ab);
            
            % 对于A->C
            obj.z_ac = alpha * obj.z_ac + (1-alpha) * 0.5 * (obj.tSA(3) - obj.tSC(2) + obj.mu_ac);
            
            % 对于B->C
            obj.z_bc = alpha * obj.z_bc + (1-alpha) * 0.5 * (obj.tSB(3) - obj.tSC(3) + obj.mu_bc);
        end

        function UpdateInnerDualVariables_Corrected(obj)
            % 修正的内层对偶变量更新
            % 使用过松弛加速
            
            omega = 1.5;  % 过松弛参数
            
            obj.mu_ab = obj.mu_ab + omega * (obj.tSA(2) - obj.tSB(2) - 2*obj.z_ab);
            obj.mu_ac = obj.mu_ac + omega * (obj.tSA(3) - obj.tSC(2) - 2*obj.z_ac);
            obj.mu_bc = obj.mu_bc + omega * (obj.tSB(3) - obj.tSC(3) - 2*obj.z_bc);
        end

        function UpdateOuterConsensus_Corrected(obj, t)
            % 修正的外层共识变量更新
            
            % 计算实际值
            S_actual = obj.ComputeActualSupply(t);
            D_actual = [sum(obj.tDA), sum(obj.tDB), obj.tDC];
            
            % 更新共识变量（使用平滑更新）
            alpha = 0.7;
            obj.s_consensus = alpha * obj.s_consensus + (1-alpha) * S_actual;
            obj.d_consensus = alpha * obj.d_consensus + (1-alpha) * D_actual;
        end

        function UpdateOuterDualVariables_Corrected(obj, t)
            % 修正的外层对偶变量更新
            
            % 计算实际值
            S_actual = obj.ComputeActualSupply(t);
            D_actual = [sum(obj.tDA), sum(obj.tDB), obj.tDC];
            
            % 1. 更新共识约束的对偶变量
            obj.u_demand = obj.u_demand + (D_actual - obj.d_consensus);
            obj.u_supply = obj.u_supply + (S_actual - obj.s_consensus);
            
            % 2. 更新价格（功率平衡约束的对偶变量）
            % 使用自适应步长
            power_imbalance = S_actual - D_actual;
            step_size = 0.1 / (1 + 0.01*obj.outer_iter_count);  % 递减步长
            
            obj.lambda = obj.lambda + step_size * power_imbalance;
            obj.lambda = max(obj.lambda, 0.01);  % 保持价格为正
        end

        function [r_norm, s_norm] = ComputeInnerResiduals_Corrected(obj, z_old)
            % 计算内层残差
            r_ab = obj.tSA(2) - obj.tSB(2) - 2*obj.z_ab;
            r_ac = obj.tSA(3) - obj.tSC(2) - 2*obj.z_ac;
            r_bc = obj.tSB(3) - obj.tSC(3) - 2*obj.z_bc;
            r_norm = norm([r_ab; r_ac; r_bc]);
            
            z_new = [obj.z_ab; obj.z_ac; obj.z_bc];
            s_norm = obj.rho_inner * norm(z_new - z_old);
        end

        function [r_norm, s_norm, power_balance] = ComputeOuterResiduals(obj, t)
            % 计算外层残差
            S_actual = obj.ComputeActualSupply(t);
            D_actual = [sum(obj.tDA), sum(obj.tDB), obj.tDC];
            
            % 原始残差：共识约束违反
            r_demand = D_actual - obj.d_consensus;
            r_supply = S_actual - obj.s_consensus;
            r_norm = norm([r_demand, r_supply]);
            
            % 对偶残差：对偶变量的变化（简化处理）
            s_norm = obj.rho_outer * (norm(obj.u_demand) + norm(obj.u_supply));
            
            % 功率平衡
            power_balance = S_actual - D_actual;
        end

        function converged = CheckOuterConvergence_Corrected(obj, r_norm, s_norm, power_balance)
            % 修正的外层收敛检查
            power_imbalance_norm = norm(power_balance);
            
            % 三个收敛条件
            converged = (r_norm < obj.eps_outer) && ...
                       (s_norm < obj.eps_outer) && ...
                       (power_imbalance_norm < obj.eps_outer);
        end

        function AdaptParameters_Corrected(obj, outer_iter, r_norm, s_norm)
            % 自适应调整参数
            
            % 每10次迭代调整一次
            if mod(outer_iter, 10) == 0
                % 基于残差比例调整
                if r_norm > 2 * s_norm && r_norm > 1e-2
                    % 原始残差太大，增加惩罚
                    obj.rho_outer = min(obj.rho_outer * 1.2, 10);
                    fprintf('  Increased ρ_outer to %.3f\n', obj.rho_outer);
                elseif s_norm > 2 * r_norm && s_norm > 1e-2
                    % 对偶残差太大，减少惩罚
                    obj.rho_outer = max(obj.rho_outer / 1.2, 0.1);
                    fprintf('  Decreased ρ_outer to %.3f\n', obj.rho_outer);
                end
            end
        end

        function AdaptInnerRho_Corrected(obj, r_norm, s_norm)
            % 自适应调整内层惩罚参数
            if r_norm > 2 * s_norm && r_norm > 1e-3
                obj.rho_inner = min(obj.rho_inner * 1.1, 10);
            elseif s_norm > 2 * r_norm && s_norm > 1e-3
                obj.rho_inner = max(obj.rho_inner / 1.1, 0.1);
            end
        end

        function S = ComputeActualSupply(obj, t)
            % 计算实际供给
            S = [obj.tSA(1) + obj.MA.PV(t) - obj.tSA(2) - obj.tSA(3) + obj.tSA(4), ...
                 obj.tSB(1) + obj.MB.PV(t) + obj.tSB(2) - obj.tSB(3) + obj.tSB(4), ...
                 obj.tSC(1) + obj.MC.PV(t) + obj.tSC(2) + obj.tSC(3) + obj.tSC(4)];
        end

        function CalculateTotalCost(obj, t)
            % 计算总成本
            obj.MA.PowerGrid = obj.tSA(1);
            obj.MB.PowerGrid = obj.tSB(1);
            obj.MC.PowerGrid = obj.tSC(1);
            obj.MA.PowerEBSS = obj.tSA(4);
            obj.MB.PowerEBSS = obj.tSB(4);
            obj.MC.PowerEBSS = obj.tSC(4);

            uA = sum(obj.U_fun(obj.tDA, obj.MA.D_it(:,t)));
            uB = sum(obj.U_fun(obj.tDB, obj.MB.D_it(:,t)));
            uC = obj.U_fun(obj.tDC, obj.MC.D_it(1,t));

            obj.TotalCost = -(uA + uB + uC) ...
                + (obj.G_fun(obj.tSA(1)) + obj.G_fun(obj.tSB(1)) + obj.G_fun(obj.tSC(1)));
        end

        function DisplayResults_Corrected(obj, t)
            % 显示结果
            fprintf('\n========== 优化结果 ==========\n');
            fprintf('外层迭代次数: %d\n', obj.outer_iter_count);
            fprintf('总内层迭代次数: %d\n', obj.total_inner_iter);
            fprintf('总成本: %.4f\n', obj.TotalCost);
            
            fprintf('\n--- 价格 ---\n');
            fprintf('λ = [%.4f, %.4f, %.4f]\n', obj.lambda(1), obj.lambda(2), obj.lambda(3));
            
            fprintf('\n--- 传输功率 ---\n');
            fprintf('A->B: %.2f kW\n', obj.tSA(2));
            fprintf('A->C: %.2f kW\n', obj.tSA(3));
            fprintf('B->C: %.2f kW\n', obj.tSB(3));
            
            % 功率平衡检查
            S = obj.ComputeActualSupply(t);
            D = [sum(obj.tDA), sum(obj.tDB), obj.tDC];
            
            fprintf('\n--- 功率平衡 ---\n');
            fprintf('供给: [%.2f, %.2f, %.2f]\n', S(1), S(2), S(3));
            fprintf('需求: [%.2f, %.2f, %.2f]\n', D(1), D(2), D(3));
            fprintf('平衡误差: %.4f\n', norm(S-D));
            
            fprintf('\n--- 对偶变量 ---\n');
            fprintf('μ_ab=%.4f, μ_ac=%.4f, μ_bc=%.4f\n', obj.mu_ab, obj.mu_ac, obj.mu_bc);
            fprintf('u_demand=[%.4f, %.4f, %.4f]\n', obj.u_demand(1), obj.u_demand(2), obj.u_demand(3));
            fprintf('u_supply=[%.4f, %.4f, %.4f]\n', obj.u_supply(1), obj.u_supply(2), obj.u_supply(3));
        end

        function PlotConvergence(obj)
            % 绘制收敛曲线
            if ~isempty(obj.outer_residuals)
                figure('Position', [100, 100, 1200, 800]);
                
                subplot(2,2,1);
                plot(1:size(obj.outer_residuals,1), obj.outer_residuals(:,1), 'b-', 'LineWidth', 2);
                hold on;
                plot(1:size(obj.outer_residuals,1), obj.outer_residuals(:,2), 'r--', 'LineWidth', 2);
                xlabel('外层迭代');
                ylabel('残差');
                legend('原始残差', '对偶残差', 'Location', 'best');
                title('外层ADMM收敛历史');
                grid on;
                set(gca, 'YScale', 'log');
                
                subplot(2,2,2);
                if ~isempty(obj.prices_history)
                    plot(1:size(obj.prices_history,1), obj.prices_history(:,1), 'b-', 'LineWidth', 2);
                    hold on;
                    plot(1:size(obj.prices_history,1), obj.prices_history(:,2), 'r--', 'LineWidth', 2);
                    plot(1:size(obj.prices_history,1), obj.prices_history(:,3), 'g-.', 'LineWidth', 2);
                    xlabel('外层迭代');
                    ylabel('价格λ');
                    legend('λ₁', 'λ₂', 'λ₃', 'Location', 'best');
                    title('价格收敛历史');
                    grid on;
                end
                
                subplot(2,2,3);
                if ~isempty(obj.inner_residuals)
                    plot(1:size(obj.inner_residuals,1), obj.inner_residuals(:,1), 'b-', 'LineWidth', 2);
                    hold on;
                    plot(1:size(obj.inner_residuals,1), obj.inner_residuals(:,2), 'r--', 'LineWidth', 2);
                    xlabel('内层迭代（累计）');
                    ylabel('残差');
                    legend('原始残差', '对偶残差', 'Location', 'best');
                    title('内层ADMM残差');
                    grid on;
                    set(gca, 'YScale', 'log');
                end
                
                subplot(2,2,4);
                bar([obj.tSA(2), obj.tSA(3), obj.tSB(3)]);
                xlabel('传输线路');
                ylabel('传输功率 (kW)');
                title('最终传输功率分配');
                grid on;
                set(gca, 'XTickLabel', {'A->B', 'A->C', 'B->C'});
            end
        end
        
        function draw(obj, t)
            % 绘图函数（保持原有功能）
            figure(1);
            subplot(3,1,1);
            x = linspace(0, 2*max(obj.tDA), 100);
            plot(x, obj.U_fun(x, obj.MA.D_it(1,t)), 'b-');
            hold on;
            plot(x, obj.lambda(1)*x, 'g-');
            plot([obj.tDA(1), obj.tDA(1)], [0, max(obj.U_fun(x, obj.MA.D_it(1,t)))], 'r-', 'LineWidth', 2);
            hold off;
            title(['区域A需求优化 (λ=', num2str(obj.lambda(1), '%.3f'), ')']);
            legend('效用函数', '价格线', '最优需求');
            grid on;
            
            subplot(3,1,2);
            x = linspace(0, 2*max(obj.tSA(1), 0.1), 100);
            plot(x, obj.G_fun(x), 'b-');
            hold on;
            plot(x, obj.lambda(1)*x, 'g-');
            plot([obj.tSA(1), obj.tSA(1)], [0, max(obj.G_fun(x))], 'r-', 'LineWidth', 2);
            hold off;
            title(['区域A发电优化']);
            legend('发电成本', '价格线', '最优发电');
            grid on;
            
            subplot(3,1,3);
            x = linspace(-2*abs(obj.tSA(4)), 2*abs(obj.tSA(4)), 100);
            if isempty(x) || all(x==0)
                x = [-1, 1];
            end
            plot(x, obj.E_fun(x, obj.MA.SOC(t), -obj.MA.BC*(1-obj.MA.SOC(t))), 'b-');
            hold on;
            plot(x, obj.lambda(1)*x, 'g-');
            plot([obj.tSA(4), obj.tSA(4)], [min(obj.E_fun(x, obj.MA.SOC(t), -obj.MA.BC*(1-obj.MA.SOC(t)))), ...
                max(obj.E_fun(x, obj.MA.SOC(t), -obj.MA.BC*(1-obj.MA.SOC(t))))], 'r-', 'LineWidth', 2);
            hold off;
            title(['区域A电池优化 (Eb=', num2str(obj.tSA(4), '%.3f'), ')']);
            legend('电池成本', '价格线', '最优充放电');
            grid on;
            
            drawnow;
        end
    end

    methods (Static)
        function f = U_fun(x, b)
            beta = 100; sigma = 1;
            if isscalar(b)
                f = beta * log(1 + sigma * x) .* (x <= b) + ...
                    beta * log(1 + sigma * b) .* (x > b);
            else
                % 处理向量输入
                f = zeros(size(x));
                for i = 1:length(x)
                    if x(i) <= b(i)
                        f(i) = beta * log(1 + sigma * x(i));
                    else
                        f(i) = beta * log(1 + sigma * b(i));
                    end
                end
            end
        end

        function f = G_fun(x)
            f = 1 * x.^2 + 0.02 * x;
        end

        function f = E_fun(x, SOC, nBmax)
            alpha = (1-SOC)/SOC*exp(4);
            f = alpha*(exp(x)-1).*(x <= 0) + (0.1*x.^2 + alpha*x).*(x > 0);
        end
    end
end