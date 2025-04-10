
classdef Community < matlab.mixin.Copyable
    
    properties
        N % array of well current species populations
        R % array of well current resource concentrations
        R0 % initial resource concentrations (1xM array)
        dNdt % consumer dynamics
        dRdt % resource dynamics
        params % struct of well params
        S0 % number of species (unchanged under optimization)
        M % number of resources
    end
    
    methods
        
        function obj = Community(init_state,params,dynamics)
            
            % init_state - {N0,R0}
            % params - struct of parameters for each well
            % dynamics - {dNdt,dRdt}
            
            [obj.N,obj.R] = deal(init_state{1},init_state{2});
            [obj.dNdt,obj.dRdt] = deal(dynamics{1},dynamics{2});
            obj.R0 = obj.R;
            obj.params = params;
            obj.S0 = length(obj.N);
            obj.M = length(obj.R);
        end
        
        function [NR,par,dydt] = PrepareWell(obj)
            
            NR = [obj.N obj.R]; % Combines N and R into single array
            
            % Compresses params, records new NR and S
            par = obj.params;
            obj.S0 = length(obj.N);

            % Generates dynamical equation using compressed params
            dN = @(NR) (obj.dNdt(NR(1:obj.S0),NR(obj.S0+1:end),par));
            dR = @(NR) (obj.dRdt(NR(1:obj.S0),NR(obj.S0+1:end),par));
            dydt = @(t,NR) ([dN(NR) dR(NR)]);
        end
        
        function [Nt,Rt] = Propagate(obj,T,dt)

            % T propagation time
            % dt timestep
            
            if obj.params.ext_rel
                xthres = obj.params.ext_thres * sum(obj.N);
            else
                xthres = obj.params.ext_thres;
            end

            t_scale = linspace(0,T,T/dt);
            % Prep and prop

            [NR0,~,dydt] = obj.PrepareWell(); %compress_consumers,compress_resources)

            options = odeset('AbsTol',1e-20);
            [~,NRarray] = ode45(@(t,NR) dydt(t,NR.').',t_scale,NR0,options);

            Nt = NRarray(:,1:obj.S0);
            Rt = NRarray(:,obj.S0+1:end);

            obj.N = NRarray(end,1:obj.S0);
            obj.R = NRarray(end,obj.S0+1:end);

            obj.N(obj.N<0)=0;
            obj.R(obj.R<0)=0;

            if sum(obj.N) < xthres
                %disp('Community collapsed on propagate')
                obj.N(:) = 0;
            else
                obj.N(obj.N < xthres) = 0; % Applies ext_thres to N
            end
        end

        function obj = Passage(obj,dilfact,refresh_resources)
            
            % tm - transfer matrix (A,B) -> old well A to new well B
            % refresh_resources - if true, resources are replenished by R0

            if obj.params.ext_rel
                xthres = obj.params.ext_thres * sum(obj.N);
            else
                xthres = obj.params.ext_thres;
            end

            if sum(obj.N) < xthres
                %disp('Community collapsed on passage')
                obj.N(:) = 0;
            else
                obj.N = obj.N ./dilfact;
                obj.N(obj.N < xthres./dilfact) = 0; % Applies ext_thres to N
            end
            
            if refresh_resources
                obj.R = obj.R0;
            end
        end
        
        function [Nf,Rf] = FindSteadyState(obj,thres,dt,max_iter,cons)
            
        

            % 1 thres - error between propagations, (0<thres<1)
            % 2 dt - prop timestep
            % 3 max_iter - max # of prop (int)
            % 4 cons - # consistency thres (int)
            % 5 ext_thres
            
            if obj.params.ext_rel
                xthres = obj.params.ext_thres * sum(obj.N);
            else
                xthres = obj.params.ext_thres;
            end

            i = 0; % prop counter
            j = 0; % cons counter
            
            while j < cons && i < max_iter && sum(obj.N)>0

                if obj.params.ext_rel
                    xthres = obj.params.ext_thres * sum(obj.N);
                else
                    xthres = obj.params.ext_thres;
                end

                % initial compositional vectors
                Ncomp0 = obj.N ./sum(obj.N,2);
                Rcomp0 = obj.R ./sum(obj.R,2);
                comp0 = [Ncomp0 Rcomp0];
                
                obj.Propagate(dt,1);
                obj.N(obj.N<xthres) = 0; % remove extinct species

                % final compositional vectors
                Ncompf = obj.N ./sum(obj.N,2);
                Rcompf = obj.R ./sum(obj.R,2);
                compf = [Ncompf Rcompf];
                
                % net compositional vectors
                diff_array = abs(compf-comp0);
                diff = max(diff_array,[],"all");
                
                i = i + 1;
                
                if diff > thres
                    j = 0;
                else
                    j = j + 1;
                end
                
                if max(obj.N,[],"all") < xthres
                    %disp("Community collapsed on FSS")
                    j = cons + 1;
                end

            end
            
            if i == max_iter % prop thres reached
                Nf = NaN(size(obj.N));
                Rf = NaN(size(obj.R));
            else % outputs
                Nf = obj.N;
                Rf = obj.R;
            end
           
            if sum(obj.N) < xthres
                %disp('Community collapsed on propagate')
                obj.N(:) = 0;
            else
                obj.N(obj.N < xthres) = 0; % Applies ext_thres to N
            end
        end

        function obj = Invade(obj,inv_info,load)
            % inv_info - struct:    c = cons prefs of inv (1xM array)
            %                       m = zero growth rate (number)
            %                       g = basal growth rate (number)
            % load - percent of population which invader will comprise (between 0 and 1)

            N_inv = load * sum(obj.N);

            if sum(obj.N) == 0
                N_inv = 0.1; % In case of community collapse
            end

            obj.params.c = [obj.params.c; inv_info.c];
            obj.params.m = [obj.params.m inv_info.m];
            obj.params.g = [obj.params.g inv_info.g];

            obj.N = [obj.N N_inv];
            obj.S0 = obj.S0 + 1;

            assert(size(obj.params.c,1)==length(obj.N))
        end

        function Perturb(obj,pmax)

            load = sum(obj.N);

            for i = 1:length(obj.N)
                perturbVal = pmax - 2*pmax*rand; % perturb by (-pmax,pmax)
                if obj.N(i)>0

                    Ni0 = obj.N(i);

                    obj.N(i) = Ni0 + (load * perturbVal);

                    if obj.params.ext_rel
                        xthres = obj.params.ext_thres * sum(obj.N);
                    else
                        xthres = obj.params.ext_thres;
                    end

                    if obj.N(i) < 0
                        obj.N(i) = xthres;
                    end
                end
            end
        end

        function Simplify(obj)

            % Remove duplicates
            N_ind = find(obj.N>0);
            c_extant = obj.params.c(obj.N>0,:);
            [~,ind,~] = unique(c_extant,'rows','stable');
            dups = setdiff(1:size(c_extant,1), ind);
            obj.N(N_ind(dups)) = 0;

            % Remove extinct species
            if obj.params.ext_rel
                xthres = obj.params.ext_thres * sum(obj.N);
            else
                xthres = obj.params.ext_thres;
            end
            present = obj.N>xthres;

            obj.S0 = sum(present);
            obj.N = obj.N(present);

            obj.params.c = obj.params.c(present,:);
            obj.params.m = obj.params.m(present);
            obj.params.g = obj.params.g(present);
        end


    end
end