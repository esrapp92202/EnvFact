%% Extinction threshold analysis

%data = rawdata;
data = rawdata(rawdata.EXT_REL == 0,:);
data = data(1:100000,:);

%{
EXT_THRESarray = [10^-6 10^-5.5 10^-5 10^-4.5 10^-4 10^-3.5 10^-3 10^-2.5 10^-2];

aug = zeros(6,3);
displ= zeros(6,3);
disr = zeros(6,3);
res = zeros(6,3);

for i = 1:length(EXT_THRESarray)
    disp(i)
    [aug_bmean,aug_minCI,aug_maxCI] = CI(data(data.EXT_THRES == EXT_THRESarray(i),:),100,'AUG');
    [disp_bmean,disp_minCI,disp_maxCI] = CI(data(data.EXT_THRES == EXT_THRESarray(i),:),100,'DISP');
    [disr_bmean,disr_minCI,disr_maxCI] = CI(data(data.EXT_THRES == EXT_THRESarray(i),:),100,'DISR');
    [res_bmean,res_minCI,res_maxCI] = CI(data(data.EXT_THRES == EXT_THRESarray(i),:),100,'RES');

    aug(i,:) = [aug_bmean,aug_minCI,aug_maxCI];
    displ(i,:) = [disp_bmean,disp_minCI,disp_maxCI];
    disr(i,:) = [disr_bmean,disr_minCI,disr_maxCI];
    res(i,:) = [res_bmean,res_minCI,res_maxCI];
end
%}

figure(1)
augline = LineErrorPlot(EXT_THRESarray,aug,'AUG');
displine = LineErrorPlot(EXT_THRESarray,displ,'DISP');
disrline = LineErrorPlot(EXT_THRESarray,disr,'DISR');
resline = LineErrorPlot(EXT_THRESarray,res,'RES');




function [bmean,minCI,maxCI] = CI(data,nboot,out)
    prop = @(table) sum(table.OUTCOME == out)/length(table.COM);
    [ci,bmeans] = bootci(nboot,prop,data);
    bmean = mean(bmeans);
    minCI = bmean-ci(1);
    maxCI = ci(2)-bmean;
end

function l = LineErrorPlot(x,data,out)
    % data: [mean1 min1 max1; mean2 min2 max2; ...]
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
        fill([x flip(x)],[(data(:,1)-data(:,2)).' flip((data(:,3)+data(:,1)).')],fcolor,'EdgeColor','none')
    else
        fill([x flip(x)],[(data(:,1)-data(:,2)).' flip((data(:,3)+data(:,1)).')],fcolor,'EdgeColor','none','FaceAlpha',.5)
    end


    l = plot(x,data(:,1),'LineWidth',2,'Color',color); % mean
    plot(x,data(:,1)-data(:,2),'--','LineWidth',1,'Color',color); %ci min
    plot(x,data(:,3)+data(:,1),'--','LineWidth',1,'Color',color); %ci max

    xlim(axes1,[0 0.01])
    ylim(axes1,[0 1])

    hold(axes1,'off');
    
    %set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
    %    [1 3 5 7 9],'XTickLabel',{'$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$'},'YTickLabel',...
    %    {'0\%','20\%','40\%','60\%','80\%','100\%'},'xscale','log');

    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','XTick',...
        [0.000001 0.00001 0.0001 0.001 0.01],'YTick',[0 0.2 0.4 0.6 0.8 1], ...
        'YTickLabel',{'0\%','20\%','40\%','60\%','80\%','100\%'},'xscale','log');

    %title('Richness-Resistance Curve','FontSize',30,...
    %    'FontName','EB Garamond');
    xlabel({'Extinction Threshold'},'FontSize',24,...
        'FontName','EB Garamond');
    ylabel({'Invasion Outcomes'},'FontSize',24,'FontName','EB Garamond');
   
end