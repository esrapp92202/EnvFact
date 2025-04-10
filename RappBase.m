% Base functions for community generation and invasion

classdef RappBase
    methods(Static)

        function iComm = InitializeCommunity(assumptions)

            import MarsBase.*

            SK = assumptions.SK;
            ENV = assumptions.ENV;
            S0 = assumptions.S0;

            Rin = zeros(1,18); % Initialize resource vector
            Rsin = randsample(18,ENV); % Determine input resource indicies
            Rin_dist = rand(1,ENV); % Set distribution of input resources
            Rin(Rsin) = Rin_dist * SK/sum(Rin_dist); % Finalize resource vector
            assumptions.R0 = Rin; % Set resource vector

            N0 = ones(1,S0);
            R0 = Rin;
            params = MarsBase.MakeParams(assumptions);
            dNdt = MakeConsumerDynamics();
            dRdt = MakeResourceDynamics();

            init_state = {N0,R0};
            dynamics = {dNdt dRdt};

            iComm = Community(init_state,params,dynamics);
            iComm.N(iComm.N < sum(iComm.N) * iComm.params.ext_thres) = 0; % initial species threshold
            iComm.Simplify;
        end

        function Comm = GenerateCommunity(assumptions)

            import RappBase.InitializeCommunity

            viable = false;
            while viable == false
                Comm = InitializeCommunity(assumptions);

                % Prop, Pass, Perturb
                Comm.Propagate(10,0.1);
                Comm.FindSteadyState(0.001,10,100,5);

                for i = 1:9
                    Comm.Passage(10,false);
                    Comm.Propagate(10,0.1);
                end

                for j = 1:10
                    if sum(Comm.N)>0
                        Comm.Perturb(0.1);
                        Comm.Propagate(10,0.1);
                        Comm.FindSteadyState(0.001,10,100,5);
                    end
                end

                if sum(Comm.N) > 0.1
                    viable = true;
                end
            end
            Comm.Simplify;
        end

        function InvadeCommunity(Comm,PropSize)

            import RappBase.GeneralConsumer

            c = GeneralConsumer(Comm.M,10,10); % Generate invader consumption profile

            while ismember(c,Comm.params.c) % Pick unique invader profile
                c = GeneralConsumer(Comm.M,10,10);
            end

            InvaderInfo = struct('M',Comm.M,'c',c,'m',1,'g',1); % Define invader
            Comm.Invade(InvaderInfo,PropSize); % Place invader

            % Find robust steady state
            Comm.Propagate(10,0.1);
            Comm.FindSteadyState(0.001,10,100,5);

            for i = 1:9
                Comm.Passage(10,false);
                Comm.Propagate(10,0.1);
            end

            for j = 1:5
                if sum(Comm.N)>0
                    Comm.Perturb(0.1);
                    Comm.Propagate(10,0.1);
                    Comm.FindSteadyState(0.001,10,100,5);
                end
            end
        end

        function c = GeneralConsumer(M,sumc,sigc)
            c_mean = sumc/M; % Average consumption per resource
            c_var = sigc^2/M; % Variance of consumption per resource
            thetac = c_var/c_mean; % Gamma distribution variable
            kc = c_mean^2/c_var; % Gamma distribution variable
            c = gamrnd(kc,thetac,1,M); % Generate c from gamma distributon
            c = c./(sum(c)/sumc); % Normalize c

            c(c<(sumc*0.2)) = 0; % Remove consumption under 20%
            c = c * sumc/sum(c); % Renormalize c
        end

    end
end