%% EXP02: Generate communities, remove one species (percolation)

assumptions0 = struct('M', 18, ...Number of resources
    'muc', 10, ...Sum of consumption rates
    'sigc', 10, ...Standard deviation used for generation of consumption rates
    'regulation','independent', ...metabolic regulation (see dRdt)
    'response','typeI', ...functional response (see dRdt)
    'supply','external', ...resource supply (see dRdt)
    'tau',ones(1,18), ...
    'w',ones(1,18),...
    'sparsity',1,...
    'l',0.8*ones(1,18)...
    );


N=100000;
exp_info = table('Size',[N 7], ...
    'VariableTypes',{'double','double','double','double','double','double','double'},...
    'VariableNames',{'COM',   'SK',    'ENV',   'RICHi', 'LOADi', 'RICHf', 'LOADf'});


c = parcluster('Processes');
c.NumWorkers = 2;
saveProfile(c);
addAttachedFiles(gcp,'RappBase.m')
%}

parfor i = 1:N

    import RappBase.*

    %% Assumption Parameters
    assumptions = assumptions0;
    ENV =      10^round(rand); % 1 or 10
    SK =       20;
    S = 50; % [5,50]

    assumptions.ENV = ENV;
    assumptions.SK = SK;

    Comm = RappBase.GenerateCommunity(assumptions,S,ENV);

    disp("Community "+string(i)+" generated")

    results = RemovalAnalysis(Comm);

    exp_info(i,:) = {i,SK,ENV,results(1),results(2),results(3),results(4)};

end

save('100k_exp02.mat', 'exp_info', '-v7.3')

%% INVASION

function results = RemovalAnalysis(Comm)

    import RappBase.*

    LOADi = sum(Comm.N);
    RICHi = sum(Comm.N>0);

    % Remove random species
    target = randi(length(Comm.N));
    Comm.N(target) = 0;
    Comm.Simplify;

    % Find robust steady state
    Comm.Propagate(10,0.1);
    Comm.FindSteadyState(0.001,5,100,5);

    for i = 1:4
        Comm.Passage(10,false);
        Comm.Propagate(10,0.1);
    end

    for j = 1:5
        if sum(Comm.N)>0
            Comm.Perturb(0.1);
            Comm.Propagate(10,0.1);
            Comm.FindSteadyState(0.001,5,100,5);
        end
    end
%}
    LOADf = sum(Comm.N);
    RICHf = sum(Comm.N>0);

    results = [RICHi, LOADi, RICHf, LOADf];

end



