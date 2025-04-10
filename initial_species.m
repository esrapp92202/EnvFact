%% Initial Species Count Analysis

% Load 2M_exp01 as rawdata

%{
S_res = zeros(50,5);
for i = 5:50

    data = rawdata(rawdata.S == i,:);
    S_res(i,:) = InitSAnalysis(data);
    
end
%}

data = rawdata(1:100000,:);
Sarray = 5:50;

aug = zeros(6,3);
disp = zeros(6,3);
disr = zeros(6,3);
res = zeros(6,3);


for i = 1:length(Sarray)
    [aug_bmean,aug_minCI,aug_maxCI] = CI(data(data.S == Sarray(i),:),100,'AUG');
    [disp_bmean,disp_minCI,disp_maxCI] = CI(data(data.S == Sarray(i),:),100,'DISP');
    [disr_bmean,disr_minCI,disr_maxCI] = CI(data(data.S == Sarray(i),:),100,'DISR');
    [res_bmean,res_minCI,res_maxCI] = CI(data(data.S == Sarray(i),:),100,'RES');

    aug(i,:) = [aug_bmean,aug_minCI,aug_maxCI];
    disp(i,:) = [disp_bmean,disp_minCI,disp_maxCI];
    disr(i,:) = [disr_bmean,disr_minCI,disr_maxCI];
    res(i,:) = [res_bmean,res_minCI,res_maxCI];
end


rich = zeros(6,3);
for i = 1:length(Sarray)
    [rich_bmean,rich_minCI,rich_maxCI] = CI(data(data.S == Sarray(i),:),100,'RICH');
    rich(i,:) = [rich_bmean,rich_minCI,rich_maxCI];
end
%}

figure(1)
augline = LineErrorPlot(aug,'AUG');
displine = LineErrorPlot(disp,'DISP');
disrline = LineErrorPlot(disr,'DISR');
resline = LineErrorPlot(res,'RES');

figure(2)
richline = LineErrorPlot(rich,'');
ylabel({'Richness'},'FontSize',24,'FontName','EB Garamond');
set(gca,'YTick',[1 2 3 4 5 6 7 8],'YTickLabel',{'1','2','3','4','5','6','7','8'},'ylim',[1 8]);



function [bmean,minCI,maxCI] = CI(data,nboot,out)
    prop = @(table) sum(table.OUTCOME == out)/length(table.COM);
    if strcmp(out,'RICH')
        prop = @(table) sum(table.RICH)/length(table.COM);
    end
    [ci,bmeans] = bootci(nboot,prop,data);
    bmean = mean(bmeans);
    minCI = bmean-ci(1);
    maxCI = ci(2)-bmean;
end

function summary_data = RichnessAnalysis(data)

    rich_sort = struct();
    summary_data = table('Size',[max(data.RICH)-min(data.RICH)+1 7], ...
        'VariableTypes',{'double','double','double','double','double','double','double'},...
        'VariableNames',{'N',     'LOAD',  'ENV',   'S',     'LOADf', 'RICHf', 'RES'});


    for i = min(data.RICH):max(data.RICH)
        rich_sort.('R'+string(i)) = data(data.RICH == i,:);

        summary_data(i-min(data.RICH)+1,:) = {length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).LOAD)/length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).ENV)/length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).S)/length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).LOADf)/length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).RICHf)/length(rich_sort.('R'+string(i)).COM),...
            sum(rich_sort.('R'+string(i)).OUTCOME == 'RES') / length(rich_sort.('R'+string(i)).COM)};
    end
end

function summary_data = InitSAnalysis(data)

    aug = sum(data.OUTCOME == 'AUG')/length(data.COM);
    disp = sum(data.OUTCOME == 'DISP')/length(data.COM);
    disr = sum(data.OUTCOME == 'DISR')/length(data.COM);
    res = sum(data.OUTCOME == 'RES')/length(data.COM);
    rich = sum(data.RICH)/length(data.COM);
    summary_data = [aug disp disr res rich];

end

function InitSpeciesPlot(data_LS,data_MS,data_HS)

    figure;

    % Create axes
    axes1 = axes;
    hold(axes1,'on');

    Y1 = data_LS;
    Y2 = data_MS;
    Y3 = data_HS;

    X1 = find(~isnan(Y1),1):find(~isnan(Y1), 1, 'last');
    X2 = find(~isnan(Y2),1):find(~isnan(Y2), 1, 'last');
    X3 = find(~isnan(Y3),1):find(~isnan(Y3), 1, 'last');

    Y1 = data_LS(~isnan(data_LS));
    Y2 = data_MS(~isnan(data_MS));
    Y3 = data_HS(~isnan(data_HS));

    % Create plot
    plot(X1,Y1,'DisplayName','Low','LineWidth',4,'Color',[0.44,0.88,0.11]);
    plot(X2,Y2,'DisplayName','Mid','LineWidth',4,'Color',[0.47,0.67,0.19]);
    plot(X3,Y3,'DisplayName','High','LineWidth',4,'Color',[0.00,0.50,0.00]);

    % Create ylabel
    ylabel({'Resistance Outcomes'},'FontSize',24,'FontName','EB Garamond');

    % Create xlabel
    xlabel({'Initial Community Richness'},'FontSize',24,...
        'FontName','EB Garamond');

    % Create title
    title({'Richness-Resistance Curve','by Initial Species Count'},'FontSize',30,...
        'FontName','EB Garamond');

    xlim(axes1,[1 13]);
    ylim(axes1,[0 1]);

    hold(axes1,'off');
    
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
        [1 3 5 7 9 11 13],'XTickLabel',...
        {'1','3','5','7','9','11','13'},'YTick',[0 0.2 0.4 0.6 0.8 1], ...
        'YTickLabel',{'0\%','20\%','40\%','60\%','80\%','100\%'});
   
    legend1 = legend(axes1,'show');
    set(legend1,'Interpreter','latex','FontSize',14);
    title(legend1,'Species');

end

function l = LineErrorPlot(data,out)
    % data: [mean1 min1 max1; mean2 min2 max2; ...]
    n = size(data,1);
    axes1 = gca;
    hold(axes1,'on');

    if strcmp(out,'AUG')
        color = [0 0.4470 0.7410];
        fcolor = color;
    elseif strcmp(out,'DISP')
        color = [0.8500 0.3250 0.0980];
        fcolor = color;
    elseif strcmp(out,'DISR')
        color = [0.4940 0.1840 0.5560];
        fcolor = color;
    elseif strcmp(out,'RES')
        color = [0.9290 0.6940 0.1250];
        fcolor = color;
    else
        color = 'k';
        fcolor = [0.8 0.8 0.8];
    end

    if isequal(fcolor,[0.8 0.8 0.8])
        fill([1:n flip(1:n)],[(data(:,1)-data(:,2)).' flip((data(:,3)+data(:,1)).')],fcolor,'EdgeColor','none')
    else
        fill([1:n flip(1:n)],[(data(:,1)-data(:,2)).' flip((data(:,3)+data(:,1)).')],fcolor,'EdgeColor','none','FaceAlpha',.5)
    end


    l = plot(1:n,data(:,1),'LineWidth',2,'Color',color); % mean
    plot(1:n,data(:,1)-data(:,2),'--','LineWidth',1,'Color',color); %ci min
    plot(1:n,data(:,3)+data(:,1),'--','LineWidth',1,'Color',color); %ci max

    xlim(axes1,[1 n])
    ylim(axes1,[0 1])

    hold(axes1,'off');
    
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
        [1 6 11 16 21 26 31 36 41 46],'XTickLabel',...
        {'5','10','15','20','25','30','35','40','45','50'},'YTickLabel',...
        {'0\%','20\%','40\%','60\%','80\%','100\%'});

    %title('Richness-Resistance Curve','FontSize',30,...
    %    'FontName','EB Garamond');
    xlabel({'Initial Species Count'},'FontSize',24,...
        'FontName','EB Garamond');
    ylabel({'Invasion Outcomes'},'FontSize',24,'FontName','EB Garamond');
   
end