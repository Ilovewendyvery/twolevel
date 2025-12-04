function visualize_results_impl(X, G, B, PV, SOC, D, time_price, ...
    convergence_info, total_cost, lambda)
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
            % 使用原电池成本函数
            B_val = B_test(j);
            SOC_val = SOC_test(i);
            if B_val <= 0
                cost = (1-SOC_val)/SOC_val * exp(0) * 0.1 * (exp(B_val) - 1);
            else
                cost = 0.00011*B_val^2 + (1-SOC_val)/SOC_val * exp(0) * 0.1 * B_val;
            end
            cost_matrix(i, j) = cost;
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
    
    % 8. 对偶变量（价格信号）
    subplot(3, 3, 8);   
    plot(1:T, lambda, 'b-', 'LineWidth', 2); hold on; 
    xlabel('时间 (小时)'); 
    ylabel('对偶变量 (价格信号)');
    legend('价格', 'Location', 'best');
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