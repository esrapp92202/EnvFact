%% EXP03: Generate communities, invade with different propagule size

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


N=50000;
propagule_sizes = logspace(-6,0,20);
exp_info = table('Size',[N*length(propagule_sizes) 8], ...
    'VariableTypes',{'double','double','double','categorical','double','double','double',  'cell'},...
    'VariableNames',{'COM',   'SK',    'ENV',   'OUTCOME',    'RICHi',  'RICHf','PROPSIZE','EXTRANK'});

c = parcluster('Processes');
c.NumWorkers = 40;
saveProfile(c);
addAttachedFiles(gcp,'RappBase.m')

parfor i = 1:N*length(propagule_sizes)

    import RappBase.*

    ENVarray = [1 2 3 5 7 10];

    %% Assumption Parameters
    assumptions = assumptions0;
    ENV =      ENVarray(randi(length(ENVarray))); % 1,2,3,5,7,10
    SK =       3*randi([1,6]); % 3 or 6 or 9 or 12 or 15 or 18
    S = 50;

    assumptions.ENV = ENV;
    assumptions.SK = SK;

    PropSize = 0.1; % propagule_sizes(ceil(i/N));
    Comm = RappBase.GenerateCommunity(assumptions,S,ENV);
    RICHi = sum(Comm.N>0);

    disp("Community "+string(i)+" generated")

    results = PropaguleAnalysis(Comm,PropSize);

    exp_info(i,:) = {i,SK,ENV,results.OUTCOME,RICHi,results.RICHf,PropSize,num2cell(results.EXTRANK,2)}; % COM,SK,ENV,OUTCOME,RICHi,RICHf

end

save('20x50k_exp03.mat', 'exp_info', '-v7.3')

%% INVASION

function results = PropaguleAnalysis(Comm,PropSize)

    import RappBase.*

    Ni = Comm.N;
    RICHi = sum(Ni>0);
    [~,poprank] = sort(Ni);

    InvadeCommunity(Comm,PropSize);

    RICHf = sum(Comm.N>0);
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

    results = struct('OUTCOME',OUTCOME,'RICHf',RICHf,'EXTRANK',EXTRANK);

end