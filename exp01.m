%% EXP01: Generate and invade communities

assumptions0 = struct('M', 18, ...Number of resources
    'muc', 10, ...Sum of consumption rates
    'sigc', 10, ...Standard deviation used for generation of consumption rates
    'regulation','independent', ...metabolic regulation (see dRdt)
    'response','typeI', ...functional response (see dRdt)
    'supply','external', ...resource supply (see dRdt)
    'tau',ones(1,18), ...
    'w',ones(1,18),...energy content per resource
    'l',0.8*ones(1,18),...leakage fraction per resource
    'spread',1, ...metabolic byproduct spread
    'Dthres',0.05, ...minimum nonzero value of metabolic matrix
    'ext_rel',false, ...absolute or relative extinction threshold
    'ext_thres',0.001 ...extinction threshold);


N = 2000000; % number of communities/invasions
exp_info = table('Size',[N 12], ...
    'VariableTypes',{'double','double','double','double','double','double','double','double',  'categorical', 'double','double','cell'},...
    'VariableNames',{'COM',   'LOAD',  'ENV',   'S0',     'RICH',  'SK',    'l',    'Spread',  'OUTCOME',     'LOADf', 'RICHf', 'EXTRANK' });

c = parcluster('Processes');
c.NumWorkers = 40;
saveProfile(c);
addAttachedFiles(gcp,'RappBase.m')

for i = 1:N

    import RappBase.*

    %% Assumption Parameters

    ENVarray = [1 2 3 5 7 10]; % Diversity domain
    SKarray = [3 6 9 12 15 18]; % Supply domain

    assumptions = assumptions0;

    % COMMUNITY VARIABLES
    ENV =      ENVarray(randi(length(ENVarray))); % Pick diversity
    SK =       SKarray(randi(length(SKarray))); % Pick supply
    S0 = randi([5,50]); % Pick S0 from 5 to 50 species

    assumptions.ENV = ENV;
    assumptions.SK = SK;
    assumptions.S0 = S0;

    Comm = RappBase.GenerateCommunity(assumptions);

    disp("Community "+string(i)+" generated")

    LOAD = sum(Comm.N);
    RICH = sum(Comm.N>0);

    results = InvasionAnalysis(Comm);
    exp_info(i,:) = {i,LOAD,ENV,assumptions.S0,RICH,SK,assumptions.l(1),assumptions.spread,results.OUTCOME,results.LOADf,results.RICHf,num2cell(results.EXTRANK,2)};

end

save('2M_exp01_S5-50', 'exp_info', '-v7.3')

%% INVASION

function results = InvasionAnalysis(Comm)

    import RappBase.*

    Ni = Comm.N;
    RICHi = sum(Ni>0);
    [~,poprank] = sort(Ni);

    InvadeCommunity(Comm,0.01);

    RICHf = sum(Comm.N>0);
    LOADf = sum(Comm.N);
    EXTRANK = RICHi + 1 - find(Comm.N(poprank) == 0);

    if isempty(EXTRANK)
        EXTRANK = 0;
    end

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

    results = struct('OUTCOME',OUTCOME,'LOADf',LOADf,'RICHf',RICHf,'EXTRANK',EXTRANK);

end



