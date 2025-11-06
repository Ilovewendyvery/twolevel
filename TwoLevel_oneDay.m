classdef TwoLevel_oneDay < handle

    properties
        PowerDemandDay;  %[D1,D2,D3]*24
        PowerSupplyDay;  %[S1,S2,S3]*24
        LambdaDay;       %[lambda1,lambda2,lambda3]*24

        %Adjustable parameters
        NT=10;
        oneHour;
        PVa;
        PVb;
        PVc;
        PBDay;
        SOCDay;

        %show

        %Other intermediate variables
        TotalCost;
        PowerGrid;mu;
        Powertrana2b;
        Powertrana2c;
        Powertranb2c;
    end

    methods
        function obj = TwoLevel_oneDay()
            obj.PowerDemandDay=zeros(obj.NT,3);
            obj.PowerSupplyDay=zeros(obj.NT,3);
            obj.LambdaDay=zeros(obj.NT,3);

            xT=linspace(0,2*pi,obj.NT);
            obj.PVa=1*(sin(xT)+1);
            obj.PVb=1*(cos(xT)+1);
            obj.PVc=1*(cos(xT)+1);

            obj.oneHour = TwoLevel_onehour();
            obj.SOCDay=zeros(obj.NT+1,3);obj.SOCDay(1,:)=zeros(1,3)+0.5;
            Solve_OneDay(obj);
        end

        function Solve_OneDay(obj)
            for k=1:obj.NT
                %Update photovoltaic power generation capacity
                obj.oneHour.PVa=obj.PVa(k);
                obj.oneHour.PVb=obj.PVb(k);
                obj.oneHour.PVc=obj.PVc(k);

                obj.oneHour.Solve();

                for kk=1:3
                    if obj.oneHour.Lambda(kk)>5
                        obj.oneHour.PB(kk)=-0.25;
                    elseif obj.oneHour.Lambda(kk)<4.5
                        obj.oneHour.PB(kk)=0.25;
                    end
                end
                obj.oneHour.Solve();


                obj.PowerDemandDay(k,:)=obj.oneHour.PowerDemand;
                obj.PowerSupplyDay(k,:)=obj.oneHour.PowerSupply;
                obj.LambdaDay(k,:)=obj.oneHour.Lambda;
                obj.Powertrana2b(k,:)=obj.oneHour.Powertrana2b;
                obj.Powertrana2c(k,:)=obj.oneHour.Powertrana2c;
                obj.Powertranb2c(k,:)=obj.oneHour.Powertranb2c;
                obj.TotalCost(k,:)=obj.oneHour.TotalCost;

                obj.PowerGrid(k,:)=obj.oneHour.PowerGrid;
                obj.PBDay(k,:)=obj.oneHour.PB;

                BC=[obj.oneHour.MA.battery_capacity,obj.oneHour.MB.battery_capacity,obj.oneHour.MC.battery_capacity];
                dt=24/obj.NT;


                obj.SOCDay(k+1,:)=obj.SOCDay(k,:)+obj.PBDay(k,:)./BC*dt;
            end
        end
        function showCombined(obj)
            % 显示六个子图的综合功率传输情况
            figure('Position', [100, 100, 1600, 1000], 'Color', 'white')  % 创建大图形窗口

            % 定义颜色方案 - 使用更专业的颜色
            colors = [0.2, 0.4, 0.8;    % 蓝色 - 微网a
                0.8, 0.2, 0.2;    % 红色 - 微网b
                0.2, 0.7, 0.3];   % 绿色 - 微网c

            lineStyles = {'-', '--', ':'};
            lineWidth = 2;

            % 子图1: 三个微网间功率传输
            subplot(4, 2, 1)
            x = 1:obj.NT;

            h1 = plot(x, obj.Powertrana2b(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.Powertrana2c(:,1), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.Powertranb2c(:,1), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('传输功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(a) 微网间功率传输', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a→b', '微网a→c', '微网b→c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图2: 三个微网电价变化
            subplot(4, 2, 2)
            h1 = plot(x, obj.LambdaDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.LambdaDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.LambdaDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('电价 (元/kWh)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(b) 各微网内部电价变化', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图3: 三个微网买电功率
            subplot(4, 2, 3)
            h1 = plot(x, obj.PowerGrid(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerGrid(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerGrid(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('买电功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(c) 微网单独买电功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图4: 三个微网发电功率
            subplot(4, 2, 4)
            h1 = plot(x, obj.PVa, 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PVb, 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PVc, 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('发电功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(d) 各微网PV发电功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图5: 三个微网供给功率
            subplot(4, 2, 5)
            h1 = plot(x, obj.PowerSupplyDay(:,1)-obj.PBDay(:,1), 'Color', colors(1,:), ...
                'LineStyle', lineStyles{1}, 'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerSupplyDay(:,2)-obj.PBDay(:,2), 'Color', colors(2,:), ...
                'LineStyle', lineStyles{2}, 'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerSupplyDay(:,3)-obj.PBDay(:,3), 'Color', colors(3,:), ...
                'LineStyle', lineStyles{3}, 'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('供给功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(e) 各微网虚拟代理商供给功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图6: 需求功率
            subplot(4, 2, 6)
            h1 = plot(x, obj.PowerDemandDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerDemandDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerDemandDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('需求功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(f) 各微网用户总需求功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % 子图7: SOC状态
            subplot(4, 2, 7)
            x_soc = 1:(obj.NT+1);

            h1 = plot(x_soc, obj.SOCDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', colors(1,:));
            hold on
            h2 = plot(x_soc, obj.SOCDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 's', 'MarkerSize', 3, 'MarkerFaceColor', colors(2,:));
            h3 = plot(x_soc, obj.SOCDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', '^', 'MarkerSize', 3, 'MarkerFaceColor', colors(3,:));
            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('SOC', 'FontSize', 11, 'FontWeight', 'normal')
            title('(g) 各微网储能SOC状态', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')
            ylim([0, 1])  % SOC通常显示在0-1范围内

            % 子图8: 充放功率 - 特别优化这个图
            subplot(4, 2, 8)

            % 使用stairs绘制，但改进样式
            h1 = stairs(x, obj.PBDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth + 0.5);
            hold on
            h2 = stairs(x, obj.PBDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth + 0.5);
            h3 = stairs(x, obj.PBDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth + 0.5);

            % 添加零线参考
            yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.5);

            hold off

            xlabel('时间步', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('充放功率 (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(h) 各微网储能充放功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'微网a', '微网b', '微网c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            

            % 添加总标题
            sgtitle('多微网能源系统运行状态综合分析', 'FontSize', 18, 'FontWeight', 'bold', ...
                'Color', [0.1, 0.2, 0.4])

            % 统一调整所有子图
            allAxes = findobj(gcf, 'Type', 'axes');
            for i = 1:length(allAxes)
                set(allAxes(i), 'FontSize', 10, 'TickDir', 'out', ...
                    'Box', 'on', 'LineWidth', 0.5);
            end

            % 调整子图间距
            set(gcf, 'Color', 'white');
            drawnow;
        end



    end
end