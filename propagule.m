%% Analysis: Invader Propagule Size (Figure S4)

% import experimental data as "rawdata" (exp02)
data = rawdata;

PropSizes = unique(data.PROPSIZE);

% Bootstrap sampling
bs_results = [PropSizes zeros(length(PropSizes),3)];
for i = 1:length(PropSizes)
    [bmean,minCI,maxCI] = CI(data(data.PROPSIZE == PropSizes(i),:),100);
    bs_results(i,2:4) = [bmean,minCI,maxCI];
end

LineErrorPlot(bs_results(:,1),bs_results(:,2:4))

function [bmean,minCI,maxCI] = CI(data,nboot)
    prop = @(table) sum(table.OUTCOME == 'RES')/length(table.COM);
    [ci,bmeans] = bootci(nboot,prop,data);
    bmean = mean(bmeans);
    minCI = bmean-ci(1);
    maxCI = ci(2)-bmean;
end

function count = CountVals(bootvect,maxrich)
    count = zeros(maxrich,1);
    [C,~,ic] = unique(bootvect);
    count(C) = accumarray(ic,1);
end

function LineErrorPlot(x,data)
    % data: [mean1 min1 max1; mean2 min2 max2; ...]
    figure1 = figure;
    axes1 = axes;
    hold(axes1,'on');

    fill([x;flip(x)],[(data(:,1)-data(:,2)); flip((data(:,3)+data(:,1)))],[0.8 0.8 0.8],'EdgeColor','none')

    plot(x,data(:,1),'LineWidth',2,'Color','k'); % mean
    plot(x,data(:,1)-data(:,2),'--','LineWidth',1,'Color','k'); %ci min
    plot(x,data(:,3)+data(:,1),'--','LineWidth',1,'Color','k'); %ci max

    xlim(axes1,[0 1])
    ylim(axes1,[0 1])

    hold(axes1,'off');
    
    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
        [0.000001 0.00001 0.0001 0.001 0.01 0.1 1],'YTick',[0 0.2 0.4 0.6 0.8 1], ...
        'YTickLabel',{'0\%','20\%','40\%','60\%','80\%','100\%'},'xscale','log');

    title('Resistance by Propagule Size','FontSize',30,...
        'FontName','EB Garamond');
    xlabel({'Propagule Size relative to Community Abundance'},'FontSize',24,...
        'FontName','EB Garamond');
    ylabel({'Resistance Outcomes'},'FontSize',24,'FontName','EB Garamond');


    xline(0.001)
    annotation(figure1,'textbox',...
        [0.531770038631077,0.164592933947773,0.123011160720622,0.104276747920013],...
        'String',{'Extinction Threshold'},...
        'FontSize',16,...
        'FontName','EB Garamond',...
        'FitBoxToText','off',...
        'EdgeColor','none');
   
end
