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
            Solve_OneDay(obj);
        end

        function Solve_OneDay(obj)
            for k=1:obj.NT
                %Update photovoltaic power generation capacity
                obj.oneHour.PVa=obj.PVa(k);
                obj.oneHour.PVb=obj.PVb(k);
                obj.oneHour.PVc=obj.PVc(k);

                obj.oneHour.Solve();

                obj.PowerDemandDay(k,:)=obj.oneHour.PowerDemand;
                obj.PowerSupplyDay(k,:)=obj.oneHour.PowerSupply;
                obj.LambdaDay(k,:)=obj.oneHour.Lambda;
                obj.Powertrana2b(k,:)=obj.oneHour.Powertrana2b;
                obj.Powertrana2c(k,:)=obj.oneHour.Powertrana2c;
                obj.Powertranb2c(k,:)=obj.oneHour.Powertranb2c;
                obj.TotalCost(k,:)=obj.oneHour.TotalCost;

                obj.PowerGrid(k,:)=obj.oneHour.PowerGrid;
            end
        end

        function showCombined(obj)
            % 显示六个子图的综合功率传输情况
            figure('Position', [100, 100, 1400, 900])  % 创建大图形窗口

            % 子图1: 三个微网间功率传输
            subplot(3, 2, 1)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.Powertrana2b(:,1), 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.Powertrana2c(:,1), 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.Powertranb2c(:,1), 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('传输功率', 'FontSize', 11)
            title('(a) 微网间功率传输', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a→b', 'a→c', 'b→c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 子图2: 三个微网电价变化
            subplot(3, 2, 2)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.LambdaDay(:,1), 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.LambdaDay(:,2), 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.LambdaDay(:,3), 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('电价', 'FontSize', 11)
            title('(b) 各微网内部电价变化', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a', 'b', 'c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 子图3: 三个微网买电功率
            subplot(3, 2, 3)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.PowerGrid(:,1), 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.PowerGrid(:,2), 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.PowerGrid(:,3), 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('买电功率', 'FontSize', 11)
            title('(c) 微网单独买电功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a', 'b', 'c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 子图4: 三个微网发电功率
            subplot(3, 2, 4)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.PVa, 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.PVb, 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.PVc, 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('发电功率', 'FontSize', 11)
            title('(d) 各微网PV发电功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a', 'b', 'c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 子图5: 三个微网供给功率
            subplot(3, 2, 5)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.PowerSupplyDay(:,1), 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.PowerSupplyDay(:,2), 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.PowerSupplyDay(:,3), 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('供给功率', 'FontSize', 11)
            title('(d)各微网虚拟代理商供给功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a', 'b', 'c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 子图5: 三个微网需求功率
            subplot(3, 2, 6)
            x = 1:obj.NT;

            % 绘制曲线，设置不同的线型和颜色
            h1 = plot(x, obj.PowerDemandDay(:,1), 'b-', 'LineWidth', 1.5);
            hold on
            h2 = plot(x, obj.PowerDemandDay(:,2), 'r--', 'LineWidth', 1.5);
            h3 = plot(x, obj.PowerDemandDay(:,3), 'g:', 'LineWidth', 1.5);
            hold off

            % 图表美化
            xlabel('时间步', 'FontSize', 11)
            ylabel('需求功率', 'FontSize', 11)
            title('(d) 各微网用户总需求功率', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'a', 'b', 'c'}, 'Location', 'best', 'FontSize', 10)
            grid on

            % 添加总标题
            sgtitle('多微网能源系统运行状态综合分析', 'FontSize', 16, 'FontWeight', 'bold')

            % 调整子图间距
            set(gcf, 'Color', 'white');  % 设置背景色为白色
        end 

    end
end