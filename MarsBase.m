classdef MarsBase
    methods(Static)

        function r = drchrnd(a,n)
            % take a sample from a dirichlet distribution
            p = length(a);
            r = gamrnd(repmat(a,n,1),1,n,p);
            r = r ./ repmat(sum(r,2),1,p);
        end

        function [N0,R0] = MakeInitialState(assumptions)
            R0 = assumptions.R0;
            N0 = ones(assumptions.S0);
        end

        function [c,D] = MakeMatrices(assumptions)
            
            % Generates consumption matrix c & metabolic byproducts matrix D
            % assumptions:  'sampling' - method for sampling c

            rng("shuffle")

            import MarsBase.*
            import RappBase.*

            M = assumptions.M; % total number of resources
            S0 = assumptions.S0;

            c = zeros(S0,M); % Will grow to (SxM)

            for i = 1:S0

                ci = GeneralConsumer(M,10,10);

                while ismember(ci,c) % Pick unique invader profile
                    ci = RappBase.GeneralConsumer(M,10,10);
                end

                c(i,:) = ci;

                
            end

            % Sampling of D

            p = ones(M,1)*1/M; ...background secretion
            D = drchrnd((p * assumptions.spread).',M).';

            if isfield(assumptions,'Dthres')
                D(D<assumptions.Dthres) = 0;

                while sum(sum(D) == 0) > 0 % Make sure no column is zeroed out
                    D = drchrnd((p * assumptions.spread).',M).';
                    D(D<assumptions.Dthres) = 0;
                end

                D = D .* (1./sum(D));
            end

        end

        function params = MakeParams(assumptions)

            import MarsBase.*

            rng("shuffle")

            % Generates initial resource array, consumer and metabolic matrices
            [~,R0] = MakeInitialState(assumptions);
            [c,D] = MakeMatrices(assumptions);

            % Prepare variables
            M = assumptions.M; %total resource count
            S0 = assumptions.S0; % total species count
            ext_thres = assumptions.ext_thres;
            ext_rel = assumptions.ext_rel;

            % Creates struct of parameters for each well
            params = struct('c',c, ...sets default parameters
                'm',ones(1,S0), ...zero-growth energy uptake per species
                'w',ones(1,M), ...energy value per resource
                'D',D, ...metabolic matrix
                'g',ones(1,S0), ...growth rate/energy value per species
                'l',zeros(1,M), ...leakage fraction per resource
                'R0',R0, ...initial resource concentrations
                'tau',1, ...
                'ext_thres',ext_thres, ...
                'ext_rel',ext_rel);

            % To change any parameter below, set it in assumptions
            default_params = ['m' 'w' 'g' 'l' 'tau' 'ext_thres' 'ext_rel'];

            for i = 1:length(default_params)
                if isfield(assumptions,default_params(i))
                    params.(default_params(i)) = assumptions.(default_params(i));
                end
            end
        end

        function dRdt = MakeResourceDynamics()

            % Inputs:   N - 1xS_tot population array
            %           R - 1xM concentration array
            % Output:   function: N,R,params -> 1xM concentration array

            % Set up resource rate equation

            J_in = @(R,params) (params.c .*R .*params.w); % (SxM) c*Rw

            J_out = @(R,params) (J_in(R.*params.l,params)*(params.D.')); % (SxM) c*Rwl *D.'
            
            dRdt = @(N,R,params) ((params.R0-R)./params.tau...
                - ((J_in(R./params.w,params)).' *(N.')).'...
                + sum(J_out(R,params)./params.w.*(N.'))...
                );
        end

        function dNdt = MakeConsumerDynamics()

            % Inputs:   N - 1xS_tot population array
            %           R - 1xM concentration array
            % Output:   function: N,R,params -> 1xS_tot population array

            % Set up consumer rate equation
            
            % [S,M] Total energy consumed per species
            J_in = @(R,params) (params.c .*R .*params.w); % (SxM) c*Rw
            
            % [S,M] Total energy used for growth per species
            J_growth = @(R,params) (J_in(R.*(1-params.l),params));

            dNdt = @(N,R,params) (params.g .*N .*((sum(J_growth(R,params),2)-params.m.').')); % (1xS) gN(c*Rw)(1-l)-m

        end
    end
end













