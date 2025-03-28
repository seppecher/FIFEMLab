classdef Mesh < handle
    properties (SetAccess = private)
        domain
        h_input
        nodes
        edges
        triangles
        building_time
    end
    methods (Access = public)

        % Constructor
        function obj = Mesh(varargin)


            % PDEtoolbox check
            if ~license('test','pde_toolbox')
                error('Mesh class needs the pde_toolbox to be installed.')
            end

            % Default Mesh
            if nargin == 0
                domain = Domain;
                obj = Mesh(domain,0.1);
                return
            end
            
            % Check correct input
            if nargin == 1 || ~isnumeric(varargin{end})
                error('Wrong input argument. Mesh needs at least two inputs and last input must be positive number.');
            end

            arg1 = varargin{1};
            h = varargin{end};

            % Shortcut the domain constuction
            if ~isa(arg1,'Domain')
                domain = Domain(varargin{1:end-1});
                obj = Mesh(domain,h);
                return
            end
            
             h = varargin{2};

            if nargin == 3
                h_bound = varargin{3};
            else 
                h_bound = h;
            end


            % General constructor
            fprintf('building mesh');
            fprintf('\n');

            % Check approximate mesh size
            Nt_exp = floor(1.5*2.31*arg1.area/h^2);

            % Big mesh alert 
            if Nt_exp > 100000  
                warning(['Approximately ' num2str(1e3*floor(Nt_exp/1000)) ' triangles will be generated.']);
                fprintf('\n');
            end



            tic;
            [gm,Ibound,domain] = arg1.geometry_matrix(h_bound);
            warning('off');
            [p,e,t]=initmesh(gm,'Hmax',h,'Hgrad',1.5,'Jiggle','mean','JiggleIter',20,'MesherVersion','R2013a');
            warning('on');
            e=e';
            e = [e(:,[1 2]) Ibound(e(:,5))];
            obj.h_input = h;
            obj.domain = domain;
            obj.nodes = p';
            obj.edges = e;
            obj.triangles = t(1:3,:)';
            obj.building_time = toc;
            fprintf('done');
            fprintf('\n');

        end

        function plot(obj,c)
            % MESH.PLOT obj.PLOT plots the Mesh object.
            if nargin == 1
                c = 'k';
            end
            t = obj.triangles;
            t = [t ones(size(t,1),1)];
            obj.domain.plot;
            hold on
            h = pdemesh(obj.nodes',obj.edges',t');
            set(h, 'Color', c);
            set(h, 'LineWidth', 0.4);
            axis equal
            hold off
        end

        
       
    end
    
    methods (Access = private)
        
      
    end
end




