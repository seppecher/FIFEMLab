classdef Mesh < handle
    properties (SetAccess = private)
        domain              % Domain object from which the mesh is built

        h_input             % Requried mesh resolution
        h_bound_input       % Requried mesh resolution at boundary               
        
        nodes               % Node positions of the mesh
        triangles           % Triangles data given as Nx3 index matrix 
        boundary_edges      % Edges at the boundary of the mesh  
        centers             % Barycenters of the triangles

        I_boundary_nodes    % Indices of boundary nodes
        I_interior_nodes    % Indices of interior nodes

        Nnodes              % Number of nodes
        Ntriangles          % Number of triangles
        Nboundary_edges     % Number of boundary edges

        h_min               % Min edge size
        h_mean              % Mean edge size
        h_max               % Max edge size
        h_std               % Santard deviation of edge sizes

        q_min               % Min quality of the mesh
        q_mean              % Mean quality of the mesh

        K                   % Default stiffness matrix
        M                   % Default mass matrix
           
        building_time       % Time used to build the mesh
        memory              % Memory used to store the mesh
    end
    properties (Access = public) % Needs private
       int1
       intx
       intxx
       P1_basis
    end
    methods (Access = public)

        % Constructor
        function obj = Mesh(varargin)
            % MESH object constuctor
            %
            % A MESH object contains a triangular mesh of a 2D domain. The
            % constructor requires the Matlab pde_toolbox. The mesher
            % algorithm is the R2013a MesherVersion of the function
            % initmesh.
            %
            % mesh = Mesh(domain,h) creates a triangular mesh of target
            % maximum resolution h.
            %
            % mesh = Mesh(domain,h,h_bound) creates a mesh of target max
            % resolution h and boundary mesh resolution h_bound.
            %
            % mesh = Mesh(domain,h,'save') or mesh = Mesh(domain,h,h_bound,'save') 
            % allows to save the large meshes and load them if previously 
            % saved.
            % 
            % Documentation: doc_mesh.mlx or doc_mesh.md.

            % PDEtoolbox check
            if ~license('test','pde_toolbox')
                error('Mesh class needs the pde_toolbox to be installed.')
            end


            opt_save = false;
            switch nargin
                case 0
                    domain = Domain;
                    obj = Mesh(domain,0.1);
                    return
                case 1
                    error('Wrong input argument. Mesh needs at least two inputs and last input must be positive number.');
                case 2
                    domain = varargin{1};
                    h = varargin{2};
                    h_bound = h;
                case 3
                    domain = varargin{1};
                    h = varargin{2};
                    if ischar(varargin{3}) && strcmp(varargin{3},'save')
                        opt_save = true;
                        h_bound = h;
                    else
                        h_bound = varargin{3};
                    end
                case 4
                    domain = varargin{1};
                    h = varargin{2};
                    h_bound = varargin{3};
                    if ischar(varargin{4}) && strcmp(varargin{4},'save')
                        opt_save = true;
                    end
            end


             % Check approximate mesh size
             Nt_exp = floor(1.5*2.31*domain.area/mean([h,h_bound])^2);

            % Procedure if save option is activiated 
             if Nt_exp >= 20000 && opt_save
                 if ~exist('saved_meshes','dir')
                     disp('build saved_meshes directory')
                     mkdir('saved_meshes');
                 end
                 disp('looking for existing mesh')
                 list = dir('saved_meshes/mesh*.mat');
                 Nlist = length(list);
                 for i = 1 : Nlist
                     meshaddress = ['saved_meshes/' list(i).name];
                     load(meshaddress,'saved_h','saved_h_bound');
                     if saved_h == h &&  saved_h_bound == h_bound
                         load(meshaddress,'saved_domain');
                         if  saved_domain == domain
                             disp('mesh found')
                             disp(['load ' meshaddress]);
                             load(meshaddress,'saved_mesh');
                             obj = saved_mesh;
                             return
                         end
                     end
                 end
                 disp('mesh not found')
                 obj = Mesh(domain,h,h_bound);
                 saved_mesh = obj;
                 saved_h = h;
                 saved_h_bound = h_bound;
                 saved_domain = domain;
                 name = ['saved_meshes/mesh' num2str(Nlist+1) '.mat' ];
                 disp('save mesh')
                 save(name,'saved_mesh','saved_h','saved_h_bound','saved_domain');
                  disp('done')
                 return
             end

            % General mesh constructor
            disp('building mesh');
            tic;
            [gm,Ibound,domain2] = domain.geometry_matrix(h_bound,true);
            [p,e,t] = initmesh(gm,'Hmax',h,'Hgrad',1.5,'Jiggle','mean','JiggleIter',20,'MesherVersion','R2013a');
            e=e';
            e = [e(:,[1 2]) Ibound(e(:,5))];
            obj.h_input = h;
            obj.h_bound_input = h_bound;
            obj.domain = domain2;
            obj.nodes = p';
            obj.boundary_edges = e;
            obj.triangles = t(1:3,:)';
            disp('building mesh properties')
            mesh_prop_builder(obj);
            if obj.q_mean<0.9
                warning('poor mesh mean quality (<0.9).')
            end
            obj.building_time = ceil(100*toc)/100;
            fprintf('done');
            fprintf('\n');
        end
        
        % Display and plot
        function disp(obj)
            % DISP method of class Mesh.
            % 
            % obj.disp displays data of a Mesh object.
            disp('Mesh object')
            disp(['nodes :           ' num2str(obj.Nnodes)])
            disp(['triangles :       ' num2str(obj.Ntriangles)])
            disp(['h (input) :       ' num2str(obj.h_input)])
            disp(['h_bound (input) : ' num2str(obj.h_bound_input)])
            disp(['min edge size :   ' num2str(obj.h_min)])
            disp(['max edge size :   ' num2str(obj.h_max)])
            disp(['mean edge size :  ' num2str(obj.h_mean)])
            disp(['edge size std :   ' num2str(obj.h_std)])
            disp(['min quality :     ' num2str(obj.q_min)])
            disp(['mean quality :    ' num2str(obj.q_mean)])
            disp(['memory (kB) :     ' num2str(obj.memory)])
            disp(['Build time (s) :  ' num2str(obj.building_time)])
            disp(' ')
        end
        function plot(obj,c)
            % PLOT method of class Mesh.
            % 
            % obj.PLOT plots a Mesh object.
            if nargin == 1
                c = 'k';
            end
            t = obj.triangles;
            t = [t ones(size(t,1),1)];
            obj.domain.plot;
            hold on
            h = pdemesh(obj.nodes',obj.boundary_edges',t');
            set(h, 'Color', c);
            set(h, 'LineWidth', 0.4);
            axis equal
            hold off
        end
        
        % Function spaces
        function u = P0(obj,varargin)
            % P0 method of class Mesh.
            % 
            % This method creates P0 functions, vector fields and tensor
            % fields on the Mesh object.
            %
            % Scalar functions:
            %
            % mesh.P0(val) creates constant scalar P0 function
            % mesh.P0(expr) creates a P0 function defined by a char
            % expression containing variable names x and/or y. 
            % mesh.P0(f) creates a P0 function defined by an handle function
            % of two variable of the form f = @(x,y) ...

            Nt = obj.Ntriangles;

            % Defaut P0 function
            if nargin == 1
                u = zeros(Nt,1);
                return
            end

            order = -1;
            argin = varargin;

            % Check if the order is specified
            if nargin >= 3 && strcmp(varargin{end-1},'order')
                order = varargin{end};
                argin = varargin(1:end-2);
            end
            
            % Default tensors
            if isempty(argin)
                u = zeros(Nt,2^order);
                return
            end

            u = obj.eval_P0(argin{1:end});
            u = tensor(u,order);
        end
        function u = P1(obj,varargin)
            % P1 method of class Mesh.
            % 
            % This method creates P1 functions, vector fields and tensor
            % fields on the Mesh object.
            %
            % Scalar functions:
            %
            % mesh.P1(val) creates constant scalar P1 function
            % mesh.P1(expr) creates a P1 function defined by a char
            % expression containing variable names x and/or y. 
            % mesh.P1(f) creates a P1 function defined by an handle function
            % of two variable of the form f = @(x,y) ...

            Np = obj.Nnodes;

            % Defaut P1 function
            if nargin == 1
                u = zeros(Np,1);
                return
            end

            order = -1;
            argin = varargin;

            % Check if the order is specified
            if nargin >= 3 && strcmp(varargin{end-1},'order')
                order = varargin{end};
                argin = varargin(1:end-2);
            end
            
            % Default tensors
            if isempty(argin)
                u = zeros(Np,2^order);
                return
            end

            u = obj.eval_P1(argin{1:end});
            u = tensor(u,order);
        end
        function b = isP0(obj,u,order)
            % isP0 method of class Mesh.
            % 
            % obj.isP0(u) checks if u is formated as a P0 function, vector
            % field or tensor field.
            %
            % obj.isP0(u,ord) checks if u is formated as a P0 tensor field 
            % of order ord.
            b = true;
            if ~isa(u,'double')
                b = false;
                return
            end
            [Nl,Nc] = size(u);
            Nt = obj.Ntriangles;
            if Nl ~= Nt
                b = false;
                return
            end
            if nargin == 3 && Nc ~= 2^order
                b = false;
            end
        end
        function b = isP1(obj,u,order)
            % isP1 method of class Mesh.
            % 
            % obj.isP1(u) checks if u is formated as a P1 function, vector
            % field or tensor field.
            %
            % obj.isP1(u,ord) checks if u is formated as a P1 tensor field 
            % of order ord.
            b = true;
            if ~isa(u,'double')
                b = false;
                return
            end
            [Nl,Nc] = size(u);
            Nn = obj.Nnodes;
            if Nl ~= Nn
                b = false;
                return
            end
            if nargin == 3 && Nc ~= 2^order
                b = false;
            end
        end

        % Function plots
        function image(obj,u,ca)
            if nargin == 2
                ca = 'auto';
            end
            Nc = size(u,2);
            s = inputname(2);
            if isreal(u)
                if Nc == 1
                    obj.image_single(u);
                    clim(ca);
                else
                    N1 = floor(sqrt(Nc));
                    N2 = ceil(Nc/N1);
                    for n = 1 : Nc
                        subplot(N1,N2,n);
                        obj.image_single(u(:,n));
                        clim(ca);
                        title([s '_{' num2str(n) '})'])
                    end
                end
            else
                if Nc == 1
                    subplot 121
                    obj.image_single(real(u));
                    clim(ca);
                    title(['Re(' s ')'])
                    subplot 122
                    obj.image_single(imag(u));
                    clim(ca);
                    title(['Im(' s ')'])
                else
                    for n = 1 : Nc
                        subplot(Nc,2,2*n-1)
                        obj.image_single(real(u(:,n)));
                        clim(ca);
                        title(['Re(' s '_{' num2str(n) '})'])
                        subplot(Nc,2,2*n)
                        obj.image_single(imag(u(:,n)));
                        clim(ca)
                        title(['Im(' s '_{' num2str(n) '})'])
                    end
                end
            end
        end
        function surf(obj,u,ca)
            if nargin == 2
                ca = 'auto';
            end
            Nc = size(u,2);
            s = inputname(2);
            if isreal(u)
                if Nc == 1
                    obj.surf_single(u);
                    clim(ca);
                else
                    N1 = floor(sqrt(Nc));
                    N2 = ceil(Nc/N1);
                    for n = 1 : Nc
                        subplot(N1,N2,n);
                        obj.surf_single(u(:,n));
                        clim(ca);
                        title([s '_{' num2str(n) '})'])
                    end
                end
            else
                if Nc == 1
                    subplot 121
                    obj.surf_single(real(u));
                    clim(ca);
                    title(['Re(' s ')'])
                    subplot 122
                    obj.surf_single(imag(u));
                    clim(ca);
                    title(['Im(' s ')'])
                else
                    for n = 1 : Nc
                        subplot(Nc,2,2*n-1)
                        obj.surf_single(real(u(:,n)));
                        clim(ca);
                        title(['Re(' s '_{' num2str(n) '})'])
                        subplot(Nc,2,2*n)
                        obj.surf_single(imag(u(:,n)));
                        clim(ca)
                        title(['Im(' s '_{' num2str(n) '})'])
                    end
                end
            end
        end
        
        % Integral calculus
        function intu = integral(obj,varargin)
            if nargin == 2
                u = varargin{1};
                if obj.isP1(u)
                    n = size(u,2);
                    intu = zeros(1,n);
                    Nn = obj.Nnodes;
                    v = obj.M*ones(Nn,1);
                    for k = 1 : n
                        intu(k) = v'*u(:,k);
                    end
                    return
                end
                if ~obj.isP0(u)
                    u = obj.P0(u);
                end
            else
                u = obj.P0(varargin);
            end
            n = size(u,2);
            intu = sum(u.*(obj.int1*ones(1,n)));
        end

        % Finite element matrices
        function K = stiffness(obj,a)
            % stiffness method of class Mesh.
            % 
            % This method builds the stiffness matrix of the P1 Finite
            % Element basis:
            %
            % K_ij = int_Omega a grad(e_j).grad(e_i)
            %
            %
            % mesh.stiffness(a) where a is a P0 matrix field or could be
            % interpreted as a P0 matrix field.
            %
            % mesh.stiffness is equivalent to mesh.stiffness(a) where a is
            % the P0 identity matrix field.
            if nargin == 1
                K = obj.K;
                return
            end
            if isscalar(a)
                K = a*obj.K;
                return
            end
            if ~obj.isP0(a,2)
                a = obj.P0(a,'order',2);
            end
            
            Nt = obj.Ntriangles;
            tri = obj.triangles;
            tp = tri';         
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3); 
            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.P1_basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            a = a(rep9,:);
            
            i1 = obj.int1(rep9);
            
            agej = [a(:,1).*alphaj(:,1) + a(:,2).*alphaj(:,2) a(:,3).*alphaj(:,1) + a(:,4).*alphaj(:,2)];
            geiagej = i1.*(alphai(:,1).*agej(:,1) + alphai(:,2).*agej(:,2));
            K = sparse(I,J,geiagej);          
        end
        function M = mass(obj,c)
            % mass method of class Mesh.
            % 
            % This method builds the mass matrix of the P1 Finite
            % Element basis:
            %
            % K_ij = int_Omega c e_j.e_i
            %
            %
            % mesh.mass(c) where a is a P0 scalar function or could be
            % interpreted as a P0 scalar function.
            %
            % mesh.mass is equivalent to mesh.mass(1).
            % 
            if nargin == 1
                M = obj.M;
                return
            end
            if isscalar(c)
                M = c*obj.M;
                return
            end
            if ~obj.isP0(c,0)
                c = obj.P0(c,'order',0);
            end

            Nt = obj.Ntriangles;
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);
            
            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.P1_basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            betamat = obj.P1_basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';
            betai = beta(rep3,:);
            betaj = beta(rep33,:);
            
            c = c(rep9,:);
            
            i1 = obj.int1(rep9);
            ix = obj.intx(rep9,:);
            ixx = obj.intxx(rep9,:);
            
            V1 = alphai(:,1).*(ixx(:,1).*alphaj(:,1) + ixx(:,2).*alphaj(:,2)) +...
                alphai(:,2).*(ixx(:,3).*alphaj(:,1) + ixx(:,4).*alphaj(:,2));
            V2 = (alphai(:,1).*ix(:,1) + alphai(:,2).*ix(:,2)).*betaj +...
                (alphaj(:,1).*ix(:,1) + alphaj(:,2).*ix(:,2)).*betai;
            V3 = betai.*betaj.*i1;
            
            aeiej = c.*(V1 + V2 + V3);
            
            M = sparse(I,J,aeiej);
        end
    end
    
    methods (Access = public) % Needs private

        function mesh_prop_builder(obj)
            p = obj.nodes;
            e = obj.boundary_edges;
            t = obj.triangles;

            obj.Nnodes = size(p,1);
            obj.Ntriangles = size(t,1);
            obj.Nboundary_edges = size(e,1);

            a = p(t(:,1),:);
            b = p(t(:,2),:);
            c = p(t(:,3),:);

            obj.centers = (a + b + c)/3;
        
            obj.I_boundary_nodes = union(e(:,1),e(:,2));
            obj.I_interior_nodes = setdiff(1:size(p,1), obj.I_boundary_nodes)';
            
            q = pdetriq(p',t');
            obj.q_mean = mean(q);
            obj.q_min = min(q);

            v1 = b - a;
            v2 = c - a;
            v3 = c - b;
            
            d1 = sqrt(v1(:,1).^2 + v1(:,2).^2);
            d2 = sqrt(v2(:,1).^2 + v2(:,2).^2);
            d3 = sqrt(v3(:,1).^2 + v3(:,2).^2);
            d = [d1 ; d2 ; d3];

            [i1,ix,ixx] = obj.integrals;
            obj.int1 = i1;
            obj.intx = ix;
            obj.intxx = ixx;
            obj.P1_basis = obj.build_basis;
              
            obj.h_min = min(d);
            obj.h_max = max(d);
            obj.h_mean = mean(d);
            obj.h_std = std(d);

            obj.K = obj.stiffness(obj.P0(1,'order',2));
            obj.M = obj.mass(obj.P0(1,'order',0));

            props = properties(obj);
            total_mem = 0;
            for ii=1:length(props)
                curr_prop = obj.(props{ii});
                s = whos('curr_prop');
                total_mem = total_mem + s.bytes;
            end
            obj.memory = ceil(total_mem/1024);
        end
        function [int1,intx,intxx] = integrals(obj)
            p = obj.nodes;
            t = obj.triangles;   
            a = p(t(:,1),:);
            b = p(t(:,2),:);
            c = p(t(:,3),:);
            u = b - a;
            v = c - a;
            d = abs(u(:,1).*v(:,2) - v(:,1).*u(:,2));
            w = u + v;
            uu = [u(:,1).*u(:,1) u(:,1).*u(:,2) u(:,2).*u(:,1) u(:,2).*u(:,2)];
            vv = [v(:,1).*v(:,1) v(:,1).*v(:,2) v(:,2).*v(:,1) v(:,2).*v(:,2)];
            uv = [u(:,1).*v(:,1) u(:,1).*v(:,2) u(:,2).*v(:,1) u(:,2).*v(:,2)];
            vu = [v(:,1).*u(:,1) v(:,1).*u(:,2) v(:,2).*u(:,1) v(:,2).*u(:,2)];
            M = 2*uu + uv + vu + 2*vv;
            int1 = d/2;
            intx = [w(:,1).*d  w(:,2).*d]/6;
            intxx = [M(:,1).*d  M(:,2).*d M(:,3).*d  M(:,4).*d]/24;       
        end
        function B = build_basis(obj)
            p = obj.nodes;
            t = obj.triangles;
            u = p(t(:,2),:) - p(t(:,1),:);
            v = p(t(:,3),:) - p(t(:,1),:);
            B = zeros(obj.Ntriangles,9);
            for j = 1 : obj.Ntriangles
                M = [0 0 1 ; u(j,:) 1 ; v(j,:) 1];
                N=inv(M);
                B(j,:)=N(:)';
            end
        end
        function val = eval_P0(obj,varargin)
            N = obj.Ntriangles;
            if nargin == 2
                input = varargin{1};
                switch class(input)
                    case 'double'
                        [Nl,Nc] = size(input);
                        if Nc == N
                            input = input';
                            Nc = Nl;
                            Nl = N;
                        end
                        switch Nl
                            case 1
                                val = ones(N,1)*input;
                                return
                            case N
                                val = input;
                            otherwise
                                if Nc == 1
                                     val = ones(N,1)*input';
                                     return
                                else
                                    error('Enable to construct a P0 function: wrong size of input.')
                                end
                        end
                    case 'function_handle'
                        eps = 0.1;
                        g = obj.centers;
                        a = obj.nodes(obj.triangles(:,1),:);
                        b = obj.nodes(obj.triangles(:,2),:);
                        c = obj.nodes(obj.triangles(:,3),:);
                        a_eps = a + eps*(g - a);
                        b_eps = b + eps*(g - b);
                        c_eps = c + eps*(g - c);


                        valg = feval(input,g(:,1),g(:,2));
                        vala = feval(input,a_eps(:,1),a_eps(:,2));
                        valb = feval(input,b_eps(:,1),b_eps(:,2));
                        valc = feval(input,c_eps(:,1),c_eps(:,2));

                        w = 1/(12*(1-eps)^2);
                        val = w*(vala + valb + valc) + (1-3*w)*valg;
                        if size(val,1) == 1
                                val = ones(N,1)*val;
                        end
                    case 'char'
                        eps = 0.1;
                        g = obj.centers;
                        a = obj.nodes(obj.triangles(:,1),:);
                        b = obj.nodes(obj.triangles(:,2),:);
                        c = obj.nodes(obj.triangles(:,3),:);
                        a_eps = a + eps*(g - a);
                        b_eps = b + eps*(g - b);
                        c_eps = c + eps*(g - c);

                        x = g(:,1);
                        y = g(:,2);
                        valg = eval(input);

                        x = a_eps(:,1);
                        y = a_eps(:,2);
                        vala = eval(input);

                        x = b_eps(:,1);
                        y = b_eps(:,2);
                        valb = eval(input);

                        x = c_eps(:,1);
                        y = c_eps(:,2);
                        valc = eval(input);

                        w = 1/(12*(1-eps)^2);
                        val = w*(vala + valb + valc) + (1-3*w)*valg;
                    
                    case 'cell'
                        Ncell = length(input);
                        val = [];
                        for k = 1 : Ncell
                            val = [val obj.eval_P0(input{k})];
                        end
                    otherwise
                        error('Input type not recognized.')

                end
            else
                val = [];
                for k = 1 : nargin-1
                    val = [val obj.eval_P0(varargin{k})];
                end
            end
        end
        function val = eval_P1(obj,varargin)
            N = obj.Nnodes;
            if nargin == 2
                input = varargin{1};
                switch class(input)
                    case 'double'
                        [Nl,Nc] = size(input);
                        if Nc == N
                            input = input';
                            Nc = Nl;
                            Nl = N;
                        end
                        switch Nl
                            case 1
                                val = ones(N,1)*input;
                                return
                            case N
                                val = input;
                            otherwise
                                if Nc == 1
                                     val = ones(N,1)*input';
                                     return
                                else
                                    error('Enable to construct a P1 function: wrong size of input.')
                                end
                        end
                    case 'function_handle'
                        p = obj.nodes;
                        val = feval(input,p(:,1),p(:,2));
                        if size(val,1) == 1
                                val = ones(N,1)*val;
                        end
                    case 'char'
                        p = obj.nodes; 
                        x = p(:,1);
                        y = p(:,2);
                        val = eval(input);
                    case 'cell'
                        Ncell = length(input);
                        val = [];
                        for k = 1 : Ncell
                            val = [val obj.eval_P1(input{k})];
                        end
                    otherwise
                        error('Input type not recognized.')

                end
            else
                val = [];
                for k = 1 : nargin-1
                    val = [val obj.eval_P0(varargin{k})];
                end
            end
        end
        function surf_single(obj,u)
            N = size(u,1);
            t = obj.triangles;
            t = [t ones(length(t),1)];
            if obj.Ntriangles == N
                pdeplot(obj.nodes',[],t','xydata',u,'xystyle','flat',...
           'zdata',u,'zstyle','discontinuous','colorbar','on');
            elseif  obj.Nnodes == N
                pdeplot(obj.nodes',[],t','xydata',u,'xystyle','interp',...
           'zdata',u,'zstyle','continuous','colorbar','on');
            else
                error('Wrong data format.')
            end
            colormap jet
            axis image
        end
        function image_single(obj,u)
            N = size(u,1);
            t = obj.triangles;
            t = [t ones(length(t),1)];
            
            if obj.Ntriangles == N
                pdeplot(obj.nodes',[],t','xydata',u,'xystyle','flat','colorbar','on');
            elseif  obj.Nnodes == N
                pdeplot(obj.nodes',[],t','xydata',u,'colorbar','on');
            else
                error('Wrong data format')
            end
            colormap jet
            axis image
        end
    end
end


% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valout = tensor(valin,order)
N = size(valin,2);
n = size(valin,1);
switch order
    case -1
        valout = valin;
    case 0
        switch N
            case 1
                valout = valin;
            otherwise
                error('Unable to constrcut a scalar function from this input.')
        end
    case 1
        switch N
            case 1
                valout = valin*[1 1];
            case 2
                valout = valin;
            otherwise
                error('Unable to constrcut a vector field from this input.')
        end

    case 2
        switch N
            case 1
                valout = valin*[1 0 0 1];
            case 2
                valout = [valin(:,1) zeros(n,2) valin(:,2)];
            case 3
                valout = [valin(:,1) valin(:,3) valin(:,3) valin(:,2)];
            case 4
                valout = valin;
            otherwise
                error('Unable to constrcut a 2nd order tensor field from this input.')
        end

    case 4
        switch N
            case 1
                I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                valout = valin(:,1)*I;
            case 2
                I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                T = [1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1];
                D = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
                valout = (valin(:,1))*(I+T) + valin(:,2)*D;
            case 4
                s11 = valin(:,1);
                s12 = valin(:,2);
                s21 = valin(:,3);
                s22 = valin(:,4);
                valout = [s11.*s11 s11.*s21  s21.*s11 s21.*s21 ...
                    s11.*s12 s11.*s22  s21.*s12 s21.*s22 ...
                    s12.*s11 s12.*s21  s22.*s11 s22.*s21 ...
                    s12.*s12 s12.*s22  s22.*s12 s22.*s22 ...
                    ];
            case 8
                s11 = valin(:,1);
                s12 = valin(:,2);
                s21 = valin(:,3);
                s22 = valin(:,4);
                t11 = valin(:,5);
                t12 = valin(:,6);
                t21 = valin(:,7);
                t22 = valin(:,8);

                valout = [s11.*t11 s11.*t21  s21.*t11 s21.*t21 ...
                    s11.*t12 s11.*t22  s21.*t12 s21.*t22 ...
                    s12.*t11 s12.*t21  s22.*t11 s22.*t21 ...
                    s12.*t12 s12.*t22  s22.*t12 s22.*t22 ...
                    ];

            case 16
                valout = valin;
            otherwise
                error('Unable to constrcut a 4th order tensor field from this input.')
        end

    otherwise
        disp(order)
        error('Invalid order')

end
end