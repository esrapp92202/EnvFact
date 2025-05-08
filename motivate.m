%% Analysis: Motivate Richness-Resistance Curves (Figure 2)

% import experimental data as "rawdata" (exp01)

% Figure 2a: resistance outcomes
data = rawdata;
env_data = zeros(13,3);

for i = 1:13
    [bmean,minCI,maxCI] = CI(data(data.RICH == i,:),100,'RES');
    env_data(i,:) = [bmean,minCI,maxCI];
end

figure(1)

LineErrorPlot(env_data,'RES')

annotation(gcf,'textbox',...
        [0.7 0.8 0.2 0.1],...
        'VerticalAlignment','middle',...
        'String',{'Supply = 12','Diversity = 5'},...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',16,...
        'FitBoxToText','off');


% Figure 2b: other outcomes
data = rawdata;
aug_data = zeros(13,3);
disp_data = zeros(13,3);
disr_data = zeros(13,3);

for i = 1:13
    [aug_bmean,aug_minCI,aug_maxCI] = CI(data(data.RICH == i,:),100,'AUG');
    [disp_bmean,disp_minCI,disp_maxCI] = CI(data(data.RICH == i,:),100,'DISP');
    [disr_bmean,disr_minCI,disr_maxCI] = CI(data(data.RICH == i,:),100,'DISR');
    aug_data(i,:) = [aug_bmean,aug_minCI,aug_maxCI];
    disp_data(i,:) = [disp_bmean,disp_minCI,disp_maxCI];
    disr_data(i,:) = [disr_bmean,disr_minCI,disr_maxCI];
end

figure(2)

augline = LineErrorPlot(aug_data,'AUG');
displine = LineErrorPlot(disp_data,'DISP');
disrline = LineErrorPlot(disr_data,'DISR');

legend1 = legend([augline displine disrline],'Augment', 'Displace', 'Disrupt');
set(legend1,'Interpreter','latex','FontSize',16,'Location','northwest');
title(legend1,'Outcome');



function [bmean,minCI,maxCI] = CI(data,nboot,out)
    prop = @(table) sum(table.OUTCOME == out)/length(table.COM);
    [ci,bmeans] = bootci(nboot,prop,data);
    bmean = mean(bmeans);
    minCI = bmean-ci(1);
    maxCI = ci(2)-bmean;
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
        [1 3 5 7 9 11 13],'XTickLabel',...
        {'1','3','5','7','9','11','13'},'YTickLabel',...
        {'0\%','20\%','40\%','60\%','80\%','100\%'});

    xlabel({'Initial Community Richness'},'FontSize',24,...
        'FontName','EB Garamond');
    ylabel({'Resistance Outcomes'},'FontSize',24,'FontName','EB Garamond');
   
end
