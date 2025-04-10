%% Analysis: Crossfeeding

% import experimental data as "rawdata" (exp01)


%% GENERATE BOOTSTRAPS

%{
%% Resource limited
data = rawdata(rawdata.SK == 6,:);
mean_limited = sum(data.OUTCOME == 'RES')/length(data.COM);

ll_lsp = data(find((data.l == 0.01).*(data.SPREAD == 0.05)),:);
ll_hsp = data(find((data.l == 0.01).*(data.SPREAD == 5)),:);
hl_lsp = data(find((data.l == 0.8).*(data.SPREAD == 0.05)),:);
hl_hsp = data(find((data.l == 0.8).*(data.SPREAD == 5)),:);

[ls1,ls2,ls3] = CI(ll_lsp,100,'res'); % low leakage, narrow
[ld1,ld2,ld3] = CI(ll_hsp,100,'res'); % low leakage, spread
[hs1,hs2,hs3] = CI(hl_lsp,100,'res'); % high leakage, narrow
[hd1,hd2,hd3] = CI(hl_hsp,100,'res'); % high leakage, spread

resource_limited = [ls1,ls2,ls3;...
    ld1,ld2,ld3;...
    hs1,hs2,hs3;...
    hd1,hd2,hd3];


%% Resource abundant
data = rawdata(rawdata.SK == 18,:);
mean_abundant = sum(data.OUTCOME == 'RES')/length(data.COM);

ll_lsp = data(find((data.l == 0.01).*(data.SPREAD == 0.05)),:);
ll_hsp = data(find((data.l == 0.01).*(data.SPREAD == 5)),:);
hl_lsp = data(find((data.l == 0.8).*(data.SPREAD == 0.05)),:);
hl_hsp = data(find((data.l == 0.8).*(data.SPREAD == 5)),:);

[ls1,ls2,ls3] = CI(ll_lsp,100,'res'); % low leakage, narrow
[ld1,ld2,ld3] = CI(ll_hsp,100,'res'); % low leakage, spread
[hs1,hs2,hs3] = CI(hl_lsp,100,'res'); % high leakage, narrow
[hd1,hd2,hd3] = CI(hl_hsp,100,'res'); % high leakage, spread

resource_abundant = [ls1,ls2,ls3;...
    ld1,ld2,ld3;...
    hs1,hs2,hs3;...
    hd1,hd2,hd3];


%% Consolidated resources
data = rawdata(rawdata.ENV == 1,:);
mean_consolidated = sum(data.OUTCOME == 'RES')/length(data.COM);

ll_lsp = data(find((data.l == 0.01).*(data.SPREAD == 0.05)),:);
ll_hsp = data(find((data.l == 0.01).*(data.SPREAD == 5)),:);
hl_lsp = data(find((data.l == 0.8).*(data.SPREAD == 0.05)),:);
hl_hsp = data(find((data.l == 0.8).*(data.SPREAD == 5)),:);

[ls1,ls2,ls3] = CI(ll_lsp,100,'res'); % low leakage, narrow
[ld1,ld2,ld3] = CI(ll_hsp,100,'res'); % low leakage, spread
[hs1,hs2,hs3] = CI(hl_lsp,100,'res'); % high leakage, narrow
[hd1,hd2,hd3] = CI(hl_hsp,100,'res'); % high leakage, spread

resource_consolidated = [ls1,ls2,ls3;...
    ld1,ld2,ld3;...
    hs1,hs2,hs3;...
    hd1,hd2,hd3];


%% Diverse resources
data = rawdata(rawdata.ENV == 10,:);
mean_diverse = sum(data.OUTCOME == 'RES')/length(data.COM);

ll_lsp = data(find((data.l == 0.01).*(data.SPREAD == 0.05)),:);
ll_hsp = data(find((data.l == 0.01).*(data.SPREAD == 5)),:);
hl_lsp = data(find((data.l == 0.8).*(data.SPREAD == 0.05)),:);
hl_hsp = data(find((data.l == 0.8).*(data.SPREAD == 5)),:);

[ls1,ls2,ls3] = CI(ll_lsp,100,'res'); % low leakage, narrow
[ld1,ld2,ld3] = CI(ll_hsp,100,'res'); % low leakage, spread
[hs1,hs2,hs3] = CI(hl_lsp,100,'res'); % high leakage, narrow
[hd1,hd2,hd3] = CI(hl_hsp,100,'res'); % high leakage, spread

resource_diverse = [ls1,ls2,ls3;...
    ld1,ld2,ld3;...
    hs1,hs2,hs3;...
    hd1,hd2,hd3];

%}

%% Plots (all data)

t = tiledlayout(2,2);
set(gcf, 'Position',  [100, 100, 800, 600])

props_limited = struct('title','Resource Limited',...
    'xaxis','',...
    'yaxis','',...
    'ylims',0.25);
props_limited.xticks = {'Low Leakage, Low Sparsity','High Leakage, Low Sparsity',... 
    'Low Leakage, High Sparsity','High Leakage, Sparsity'};
DotErrorPlot(resource_limited, props_limited)

props_abundant = struct('title','Resource Abundant',...
    'xaxis','',...
    'yaxis','',...
    'ylims',0.35);
props_abundant.xticks = {'Low Leakage, Low Sparsity','High Leakage, Low Sparsity',... 
    'Low Leakage, High Sparsity','High Leakage, Sparsity'};
DotErrorPlot(resource_abundant, props_abundant)

props_consolidated = struct('title','Consolidated Resources',...
    'xaxis','',...
    'yaxis','',...
    'ylims',0.25);
props_consolidated.xticks = {'Low Leakage, Low Sparsity','High Leakage, Low Sparsity',... 
    'Low Leakage, High Sparsity','High Leakage, Sparsity'};
DotErrorPlot(resource_consolidated, props_consolidated)

props_diverse = struct('title','Diverse Resources',...
    'xaxis','',...
    'yaxis','',...
    'ylims',0.25);
props_diverse.xticks = {'Low Leakage, Low Sparsity','High Leakage, Low Sparsity',... 
    'Low Leakage, High Sparsity','High Leakage, Sparsity'};
DotErrorPlot(resource_diverse, props_diverse)

ylabel(t,'Resistance Outcomes','FontSize',24,'FontName','EB Garamond')
%title(t,'Resistance due to Consumption Leakage and Metabolic Sparsity','FontSize',24,'FontName','EB Garamond')

function [bmean,minCI,maxCI] = CI(data,nboot,type)
    prop = struct('res',@(table) sum(table.OUTCOME == 'RES')/length(table.COM), ...
        'rich',@(table) sum(table.RICH)/length(table.COM));
    [ci,bmeans] = bootci(nboot,prop.(type),data);
    bmean = mean(bmeans);
    minCI = bmean-ci(1);
    maxCI = ci(2)-bmean;
end

function g = DotErrorPlot(data,props)
    % data: [mean1 min1 max1; mean2 min2 max2; ...]
    n = size(data,1);

    %figure1 = figure;
    axes1 = nexttile;
    hold(axes1,'on');

    m = mean(data(:,1));
    
    
    xlabels = {['$$\begin{array}{c}' '{\rm{Low\;Leakage}}' '\\' '{\rm{Low\;Spread}}' '\end{array}$$'] ...
        ['$$\begin{array}{c}' '{\rm{Low\;Leakage}}' '\\' '{\rm{High\;Spread}}' '\end{array}$$'] ...
        ['$$\begin{array}{c}' '{\rm{High\;Leakage}}' '\\' '{\rm{Low\;Spread}}' '\end{array}$$'] ...
        ['$$\begin{array}{c}' '{\rm{High\;Leakage}}' '\\' '{\rm{High\;Spread}}' '\end{array}$$']};
%}


    g = errorbar(1:n,data(:,1),data(:,2),data(:,3),...
        'MarkerSize',7,'MarkerFaceColor',[0 0 0],...
        'MarkerEdgeColor','none',...
        'Marker','square',...
        'LineStyle','none',...
        'Color',[0 0 0],...
        'CapSize',18);

    xlim(axes1,[0.75 n+0.25])
    ylim(axes1,[0 1])

    title(props.title);
    ylabel(props.yaxis);
    yline(m,'--')

    set(axes1, 'FontSize', 18, 'FontName', 'EB Garamond', ...
    'XTick', 1:n, 'XTickLabel', xlabels, 'TickLabelInterpreter', 'latex', ...
    'YTick', [0 0.25 0.5 0.75 1], ...
    'YTickLabel', {'0\%', '25\%', '50\%', '75\%', '100\%'});
    axes1.XAxis.FontSize = 12;

    xticks(1:n); % Set tick positions
    xticklabels(xlabels); % Apply multiline labels
    xtickangle(axes1,40)
end



