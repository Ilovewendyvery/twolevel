classdef TwoLevel_oneDay < handle

    properties
        PowerDemandDay;  %[D1,D2,D3]*24
        PowerSupplyDay;  %[S1,S2,S3]*24
        LambdaDay;       %[lambda1,lambda2,lambda3]*24

        %Adjustable parameters
        NT=12;
        oneHour;
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
            obj.mu=zeros(obj.NT,3); 

            obj.oneHour = TwoLevel_onehour(obj.NT);
            obj.SOCDay=zeros(obj.NT+1,3);obj.SOCDay(1,:)=zeros(1,3)+0.5;
            Solve_OneDay(obj);
        end

        function Solve_OneDay(obj)
            for k=1:obj.NT 
                obj.oneHour.Solve(k); 

                obj.PowerDemandDay(k,:)=[sum(obj.oneHour.tDA), sum(obj.oneHour.tDB), sum(obj.oneHour.tDC)];
                S = [obj.oneHour.tSA(1) + obj.oneHour.MA.PV(k) - obj.oneHour.tSA(2) - obj.oneHour.tSA(3) + obj.oneHour.tSA(4), ...
                obj.oneHour.tSB(1) + obj.oneHour.MB.PV(k) + obj.oneHour.tSB(2) - obj.oneHour.tSB(3) + obj.oneHour.tSB(4), ...
                obj.oneHour.tSC(1) + obj.oneHour.MC.PV(k) + obj.oneHour.tSC(2) + obj.oneHour.tSC(3) + obj.oneHour.tSC(4)];

                obj.PowerSupplyDay(k,:)=S;
                obj.LambdaDay(k,:)=obj.oneHour.lambda;
                obj.mu(k,:)=[obj.oneHour.mu_ab; obj.oneHour.mu_ac; obj.oneHour.mu_bc];
                obj.Powertrana2b(k,:)=obj.oneHour.tSA(2);
                obj.Powertrana2c(k,:)=obj.oneHour.tSA(3);
                obj.Powertranb2c(k,:)=obj.oneHour.tSB(3);
                obj.TotalCost(k,:)=obj.oneHour.TotalCost;

                obj.PowerGrid(k,:)=[obj.oneHour.MA.PowerGrid,obj.oneHour.MB.PowerGrid,obj.oneHour.MC.PowerGrid];
                obj.PBDay(k,:)=[obj.oneHour.MA.PowerEBSS,obj.oneHour.MB.PowerEBSS,obj.oneHour.MC.PowerEBSS]; 

                BC=[obj.oneHour.MA.BC,obj.oneHour.MB.BC,obj.oneHour.MC.BC];
                dt=24/obj.NT;


                obj.SOCDay(k+1,:)=obj.SOCDay(k,:)-obj.PBDay(k,:)./BC;
                obj.oneHour.MA.SOC(k+1)=obj.SOCDay(k+1,1);
                obj.oneHour.MB.SOC(k+1)=obj.SOCDay(k+1,2);
                obj.oneHour.MC.SOC(k+1)=obj.SOCDay(k+1,3);
            end
        end
        function showCombined(obj)
            % Display comprehensive power transfer situation of six subplots
            figure('Position', [100, 100, 1600, 900], 'Color', 'white')  % Create large figure window

            % Define color scheme - using more professional colors
            colors = [0.2, 0.4, 0.8;    % Blue - Microgrid a
                0.8, 0.2, 0.2;    % Red - Microgrid b
                0.2, 0.7, 0.3];   % Green - Microgrid c

            lineStyles = {'-', '--', ':'};
            lineWidth = 2;

            % Subplot 1: Power transfer between three microgrids
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

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Transfer Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(a) Inter-microgrid Power Transfer', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a→b', 'MG a→c', 'MG b→c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 2: Electricity price changes in three microgrids
            subplot(4, 2, 2)
            h1 = plot(x, obj.LambdaDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.LambdaDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.LambdaDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Price (Yuan/kWh)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(b) Internal Electricity Price Changes', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 3: Power purchase from grid for three microgrids
            subplot(4, 2, 3)
            h1 = plot(x, obj.PowerGrid(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerGrid(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerGrid(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Grid Purchase Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(c) Individual Grid Purchase Power', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 4: Power generation in three microgrids
            subplot(4, 2, 4)
            h1 = plot(x, obj.oneHour.MA.PV, 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.oneHour.MB.PV, 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.oneHour.MC.PV, 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Generation Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(d) PV Generation Power', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 5: Supply power in three microgrids
            subplot(4, 2, 5)
            h1 = plot(x, obj.PowerSupplyDay(:,1), 'Color', colors(1,:), ...
                'LineStyle', lineStyles{1}, 'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerSupplyDay(:,2), 'Color', colors(2,:), ...
                'LineStyle', lineStyles{2}, 'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerSupplyDay(:,3), 'Color', colors(3,:), ...
                'LineStyle', lineStyles{3}, 'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Supply Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(e) Virtual Agent Supply Power', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 6: Demand power
            subplot(4, 2, 6)
            h1 = plot(x, obj.PowerDemandDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.PowerDemandDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.PowerDemandDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Demand Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(f) Total User Demand Power', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 7: SOC status
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

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('SOC', 'FontSize', 11, 'FontWeight', 'normal')
            title('(g) Energy Storage SOC Status', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')
            ylim([0, 1])  % SOC typically displayed in 0-1 range

            % Subplot 8: Charge/discharge power - specially optimized for this plot
            subplot(4, 2, 8)

            % Use stairs for plotting, but with improved style
            h1 = stairs(x, obj.PBDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth + 0.5);
            hold on
            h2 = stairs(x, obj.PBDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth + 0.5);
            h3 = stairs(x, obj.PBDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth + 0.5);

            % Add zero line reference
            yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.5);

            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Charge/Discharge Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(h) Energy Storage Charge/Discharge Power', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a', 'MG b', 'MG c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Add overall title
            sgtitle('Comprehensive Analysis of Multi-Microgrid Energy System Operation', 'FontSize', 18, 'FontWeight', 'bold', ...
                'Color', [0.1, 0.2, 0.4])

            % Uniformly adjust all subplots
            allAxes = findobj(gcf, 'Type', 'axes');
            for i = 1:length(allAxes)
                set(allAxes(i), 'FontSize', 10, 'TickDir', 'out', ...
                    'Box', 'on', 'LineWidth', 0.5);
            end

            % Adjust subplot spacing
            set(gcf, 'Color', 'white');
            drawnow;
        end

        function showCombined2(obj)
            % Display comprehensive power transfer situation of six subplots
            figure('Color', 'white')  % Create large figure window

            % Define color scheme - using more professional colors
            colors = [0.2, 0.4, 0.8;    % Blue - Microgrid a
                0.8, 0.2, 0.2;    % Red - Microgrid b
                0.2, 0.7, 0.3];   % Green - Microgrid c

            lineStyles = {'-', '--', ':'};
            lineWidth = 2;

            % Subplot 1: Power transfer between three microgrids
            subplot(2, 1, 1)
            x = 1:obj.NT;

            h1 = plot(x, obj.Powertrana2b(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.Powertrana2c(:,1), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.Powertranb2c(:,1), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Transfer Power (kW)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(a) Inter-microgrid Power Transfer', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3], {'MG a→b', 'MG a→c', 'MG b→c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off')
            %grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

            % Subplot 2: Electricity price changes in three microgrids
            subplot(2, 1, 2)
            h1 = plot(x, obj.LambdaDay(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            hold on
            h2 = plot(x, obj.LambdaDay(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth', lineWidth, 'Marker', 'none');
            h3 = plot(x, obj.LambdaDay(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', lineWidth, 'Marker', 'none');

            h4 = plot(x, obj.mu(:,1), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
                'LineWidth', 3, 'Marker', 'none');
           
            h5 = plot(x, obj.mu(:,2), 'Color', colors(2,:), 'LineStyle', lineStyles{2}, ...
                'LineWidth',3, 'Marker', 'none');
            h6 = plot(x, obj.mu(:,3), 'Color', colors(3,:), 'LineStyle', lineStyles{3}, ...
                'LineWidth', 3, 'Marker', 'none');
            hold off

            xlabel('Time Step', 'FontSize', 11, 'FontWeight', 'normal')
            ylabel('Price (Yuan/kWh)', 'FontSize', 11, 'FontWeight', 'normal')
            title('(b) Internal Electricity Price Changes', 'FontSize', 12, 'FontWeight', 'bold')
            legend([h1, h2, h3, h4, h5, h6], {'MG a', 'MG b', 'MG c','MG a→b', 'MG a→c', 'MG b→c'}, ...
                'Location', 'best', 'FontSize', 9, 'Box', 'off') 
            grid on
            set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--')

 
            % Add overall title
            sgtitle('Comprehensive Analysis of Multi-Microgrid Energy System Operation', 'FontSize', 18, 'FontWeight', 'bold', ...
                'Color', [0.1, 0.2, 0.4])

            % Uniformly adjust all subplots
            allAxes = findobj(gcf, 'Type', 'axes');
            for i = 1:length(allAxes)
                set(allAxes(i), 'FontSize', 10, 'TickDir', 'out', ...
                    'Box', 'on', 'LineWidth', 0.5);
            end

            % Adjust subplot spacing
            set(gcf, 'Color', 'white');
            drawnow;
        end



    end
end