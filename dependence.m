%% Analysis: Dependence (inverse percolation)

% import experimental data as "rawdata" (exp02)
%{
consolidated = rawdata(rawdata.ENV == 1,:);
dispersed = rawdata(rawdata.ENV == 10,:);

cons_rich_loss = consolidated.RICHi-consolidated.RICHf;
disp_rich_loss = dispersed.RICHi-dispersed.RICHf;
%}
% True distribution
%{
[cons_counts,cons_vals] = groupcounts(cons_rich_loss);
[disp_counts,disp_vals] = groupcounts(disp_rich_loss);

cons_prop = cons_counts(1:8)/sum(cons_counts);
disp_prop = disp_counts(1:8)/sum(disp_counts);
%}

% Bootstrap distribution
%{
[cons_bmean,cons_minCI,cons_maxCI] = CI(cons_rich_loss,100);
[disp_bmean,disp_minCI,disp_maxCI] = CI(disp_rich_loss,100);
%}

barcolor = [0.00 0.45 0.74; 0.30 0.75 0.93];

t = tiledlayout(1,2);
set(gcf, 'Position',  [100, 100, 1000, 500])
xlabel(t,'Total Richness Loss','FontSize',24,'FontName','EB Garamond');
ylabel(t,'Proportion of Outcomes','FontSize',24,'FontName','EB Garamond')

BarErrorPlot(cons_bmean,cons_minCI,cons_maxCI,10,barcolor(1,:));
title({'Consolidated Resources'},'FontSize',30,...
        'FontName','EB Garamond');

BarErrorPlot(disp_bmean,disp_minCI,disp_maxCI,10,barcolor(2,:));
title({'Dispersed Resources'},'FontSize',30,...
        'FontName','EB Garamond');











function [bmean,minCI,maxCI] = CI(rich_loss,nboot)

    [ci,bmeans] = bootci(nboot,@(x) CountVals(x,max(rich_loss)),rich_loss);
    bmean = mean(bmeans)/length(rich_loss);
    minCI = bmean-ci(1,:)/length(rich_loss);
    maxCI = ci(2,:)/length(rich_loss)-bmean;

end

function count = CountVals(bootvect,maxrich)
    count = zeros(maxrich,1);
    [C,~,ic] = unique(bootvect);
    count(C) = accumarray(ic,1);
end

function BarErrorPlot(mean,minCI,maxCI,richmax,barcolor)

    axes1 = nexttile;
    bar(1:richmax,mean(1:richmax),'FaceColor',barcolor)

    hold on

    errorbar(1:richmax,mean(1:richmax),minCI(1:richmax),maxCI(1:richmax),...
        'LineStyle','none','LineWidth',1,...
        'Color',[0 0 0]);

    xlim(axes1,[0 11]);
    ylim(axes1,[0 0.5]);

    set(axes1,'FontSize',20,'TickLabelInterpreter','latex','YTick',[0 0.1 0.2 0.3 0.4 0.5],'YTickLabel',...
        {'0\%','10\%','20\%','30\%','40\%','50\%'});


end