%% Analysis: Resistance by Diversity with Supply Curves (Figure 4b)

% import experimental data as "rawdata" (exp01)

data = rawdata;
res_total = struct();

for i = [3 6 9 12 15 18]
    data = rawdata(rawdata.SK == i,:);
    env_data = ResistanceAnalysis(data);
    res_total.('K'+string(i)) = env_data.RES([1 2 3 5 7 10]);
end

ENVPlot(res_total)



function ENV_summary_data = ResistanceAnalysis(data)

    ENV_sort = struct();

    ENV_summary_data = table('Size',[6 7], ...
        'VariableTypes',{'double','double','double','double','double','double','double'},...
        'VariableNames',{'N',     'LOAD',  'ENV',   'S',     'LOADf', 'RICHf', 'RES'});

    for i = [1 2 3 5 7 10]
        ENV_sort.('E'+string(i)) = data(data.ENV == i,:);

        ENV_summary_data(i,:) = {length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).LOAD)/length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).ENV)/length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).S)/length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).LOADf)/length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).RICHf)/length(ENV_sort.('E'+string(i)).COM),...
            sum(ENV_sort.('E'+string(i)).OUTCOME == 'RES') / length(ENV_sort.('E'+string(i)).COM) };
    end
end

function ENVPlot(data_struct)

    %CREATEFIGURE(Y1, Y2, Y3, YMatrix1, Y4)
    %  Y1:  vector of plot y data
    %  Y2:  vector of plot y data
    %  Y3:  vector of plot y data
    %  YMATRIX1:  matrix of plot y data
    %  Y4:  vector of plot y data

    %  Auto-generated by MATLAB on 06-Feb-2025 18:42:31

    % Create figure
    figure;

    % Create axes
    axes1 = axes;
    hold(axes1,'on');

    X = [1 2 3 5 7 10];

    Y1 = data_struct.K3;
    Y2 = data_struct.K6;
    Y3 = data_struct.K9;
    Y4 = data_struct.K12;
    Y5 = data_struct.K15;
    Y6 = data_struct.K18;

    % Create plot
    plot(X,Y1,'DisplayName','3','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.96078431372549 0.87843137254902 0.901960784313726],'Marker','.','MarkerSize',30);
    plot(X,Y2,'DisplayName','6','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.941176470588235 0.729411764705882 0.768627450980392],'Marker','.','MarkerSize',30);
    plot(X,Y3,'DisplayName','9','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.92156862745098 0.568627450980392 0.63921568627451],'Marker','.','MarkerSize',30);
    plot(X,Y4,'DisplayName','12','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.92156862745098 0.388235294117647 0.47843137254902],'Marker','.','MarkerSize',30);
    plot(X,Y5,'DisplayName','15','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.901960784313726 0.188235294117647 0.32156862745098],'Marker','.','MarkerSize',30);
    plot(X,Y6,'DisplayName','18','MarkerSize',20,'Marker','.','LineWidth',2,'Color',[0.63921568627451 0.0784313725490196 0.180392156862745],'Marker','.','MarkerSize',30);
    
    % Create ylabel
    ylabel({'Resistance Outcomes'},'FontSize',24,'FontName','EB Garamond');

    % Create xlabel
    xlabel({'# of Input Resources'},'FontSize',24,...
        'FontName','EB Garamond');

    % Create title
    %title({'Richness-Resistance Curve','by Resource Diversity'},'FontSize',30,...
    %    'FontName','EB Garamond');

    xlim(axes1,[1 10]);
    ylim(axes1,[0 1]);

    hold(axes1,'off');
    
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
        [1 2 3 5 7 10],'XTickLabel',...
        {'1','2','3','5','7','10'},'YTick',[0 0.2 0.4 0.6 0.8 1],'YTickLabel',...
        {'0\%','20\%','40\%','60\%','80\%','100\%'});
   
    legend1 = legend(axes1,'show');
    set(legend1,'Interpreter','latex','FontSize',16,'Location','southwest');
    title(legend1,'Supply');

end
