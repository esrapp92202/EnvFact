%% EXP04: Generate and invade communities with varying S

assumptions0 = struct('M', 18, ...Number of resources
    'muc', 10, ...Sum of consumption rates
    'sigc', 10, ...Standard deviation used for generation of consumption rates
    'regulation','independent', ...metabolic regulation (see dRdt)
    'response','typeI', ...functional response (see dRdt)
    'supply','external', ...resource supply (see dRdt)
    'tau',ones(1,18), ...
    'w',ones(1,18)...
    );


N=1000000;
exp_info = table('Size',[N 11], ...
    'VariableTypes',{'double','double','double','double','double','double','double','double',  'categorical', 'double','double'},...
    'VariableNames',{'COM',   'LOAD',  'ENV',   'S',     'RICH',  'SK',    'l',     'Sparsity','OUTCOME',     'LOADf', 'RICHf'});

c = parcluster('Processes');
c.NumWorkers = 40;
saveProfile(c);
addAttachedFiles(gcp,'RappBase.m')

parfor i = 1:N

    import RappBase.*

    %% Assumption Parameters

    assumptions = assumptions0;
    ENV =      5; % 1,2,3,5,7,10
    SK =       12; % 3 or 6 or 9 or 12 or 15 or 18
    l =        0.8;
    sparsity = 1;
    S = randi([5,50]); % [5,50]

    assumptions.ENV = ENV;
    assumptions.SK = SK;
    assumptions.l = l * ones(1,18);
    assumptions.sparsity = sparsity;

    Comm = RappBase.GenerateCommunity(assumptions,S,ENV);

    disp("Community "+string(i)+" generated")

    LOAD = sum(Comm.N);
    ENV = sum(Comm.R0>0);
    RICH = sum(Comm.N>0);

    results = InvasionAnalysis(Comm);
    exp_info(i,:) = {i,LOAD,ENV,S,RICH,SK,l,sparsity,results.OUTCOME,results.LOADf,results.RICHf};

end

save('1M_exp04', 'exp_info', '-v7.3')

%% INVASION

function results = InvasionAnalysis(Comm)

    import RappBase.*

    Ni = Comm.N;
    RICHi = sum(Ni>0);

    InvadeCommunity(Comm,0.01);

    RICHf = sum(Comm.N>0);
    LOADf = sum(Comm.N);

    invader_present = Comm.N(end)>0;

    if RICHi < RICHf && invader_present
        OUTCOME = 'AUG';
        disp('Outcome: AUGMENT')
    elseif RICHi >= RICHf && invader_present
        OUTCOME = 'DISP';
        disp('Outcome: DISPLACE')
    elseif RICHi > RICHf && ~invader_present
        OUTCOME = 'DISR';
        disp('Outcome: DISRUPT')
    elseif RICHi == RICHf && ~invader_present
        OUTCOME = 'RES';
        disp('Outcome: RESIST')
    end

    results = struct('OUTCOME',OUTCOME,'LOADf',LOADf,'RICHf',RICHf);

end



