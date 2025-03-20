classdef Domain

    properties %(SetAccess = private)
        nodes
        edges
        poly
        area
    end

    methods (Access = public)

        function obj = Domain(varargin)
            % DOMAIN object constuctor.
            %
            % A Domain object contains geometrical information on a 2D
            % needed to generate later a Mesh object.
            %
            % Basic domains:
            % domain = Domain or domain = Domain('square') defines the unit
            % square domain.
            %
            % General domain consctuction:
            % domain = Domain(nodes,edges) where nodes is a Nx2 matrix of N
            % nodes positions and edges is a Px1 cell defining the domain
            % boundary edges.
            

            % Default domain
            if nargin == 0
                obj = Domain([0 0 ; 1 0 ; 1 1 ; 0 1]);
                return
            end

            arg1 = varargin{1};

            % Shortcut constructors
            if ischar(arg1)
                if nargin == 1
                    obj = Domain(arg1,1);
                    return
                end
                r = varargin{2};
                switch arg1
                    case {'square','Square'}
                        obj = Domain([0 0 ; r 0 ; r r ; 0 r]);
                        return
                    case {'disc','Disc'}
                        edgefun = @(t) [r*cos(t) r*sin(t)];
                        obj = Domain([r 0],{{1,1,edgefun,[0,2*pi]}});
                        return
                    case {'image','Image'}
                        h = 5;
                        I = double(varargin{2});
                        if size(I,3) == 3
                            I = sum(I,3)/3;
                        end
                        I(I>0) = 1;
                        if nargin == 3
                            dx = varargin{3};
                        else
                            dx = 1;
                        end
                        B = bwboundaries(I);
                        [Ny,Nx] = size(I);
                        Ly = Ny*dx;
                        Nb = length(B);
                        Ind = zeros(Nb,1);
                        for n = 1 : Nb
                            b = B{n};
                            b = under_sample(b,h);
                            p = dx*[b(:,2) Ly - b(:,1)];
                            N = size(p,1);
                            B{n} = p;
                            if N>2
                                Ind(n) = 1;
                            end
                        end
                        B = B(logical(Ind));
                        Nb = length(B);

                        P = zeros(Nb,2);
                        E = cell(Nb,1);

                        for n = 1 : Nb
                            p = B{n};
                            if n == 1
                                side = 'L';
                            else
                                side = 'R';
                            end
                            P(n,:) = p(1,:);
                            E{n} = {n,n,[p ; p(1,:)],side};
                        end
                        obj = Domain(P,E);
                        return


                    otherwise
                        error('Unvalid first argument.')
                end
            end


            % Constructor from nodes only
            if nargin == 1
                nodes = arg1;
                if size(nodes,2) ~= 2
                    nodes = nodes';
                    if size(nodes,2) ~=2
                        error('Wrong format for first input: require a Nx2 nodes matrix.');
                    end
                end
                Nn = size(nodes,1);
                edges = cell(Nn,5);
                for n = 1 : Nn-1
                    edges{n,1} = n;
                    edges{n,2} = n+1;
                    edges{n,3} = 'straight';
                    edges{n,5} = 'L';
                end
                edges{Nn,1} = Nn;
                edges{Nn,2} = 1;
                edges{Nn,3} = 'straight';
                edges{Nn,5} = 'L';
                obj.nodes = nodes;
                obj.edges = edges;
                obj.poly = obj.build_poly;
                obj.area = area(obj.poly);
                obj = obj.check_structure;
                return
            end



            % Basic constructor
            obj.nodes = arg1;
            edges = varargin{2};

            [Ne,p] = size(edges);

            if p>1 && Ne == 1
                edges = edges';
            end
            [Ne,p] = size(edges);


            if p == 5
                obj.edges = edges;
                obj = obj.check_structure;
                return
            elseif p == 1
                C = cell(Ne,5);
                for n = 1 : Ne
                    C(n,:) = format_edge(edges{n});
                end
                obj.edges = C;
                obj.poly = obj.build_poly;
                obj.area = area(obj.poly);
                obj = obj.check_structure;


            else
                error('Wrong format for second input: requires a Nx1 cell.');

            end
        end

        function disp(obj)
            % DISP method: display informations contained in the Domain
            % object.
            disp('Domain object')
            disp(['area = ' num2str(obj.area)]);

            Ne = size(obj.edges,1);
            if Ne <= 50 % Full description



                disp(['nodes (' num2str(size(obj.nodes,1)) '):']);
                disp(obj.nodes)
                disp(['edges (' num2str(Ne) '):']);

                for i = 1 : Ne
                    edge = obj.edges(i,:);
                    type = edge{3};
                    side = edge{5};
                    a = num2str(edge{1});
                    b = num2str(edge{2});
                    type_expr = '';
                    switch side
                        case 'R'
                            type_expr = ', hole';
                        case 'LR'
                            type_expr = ', internal';
                    end

                    disp(['edge ' num2str(i) ': from ' a ' to ' b ' ' type type_expr]);


                end
            else % Synthetic description
                disp(['Nodes number = ' num2str(size(obj.nodes,1))]);
                disp(['Edges number = ' num2str(length(obj.edges))]);
            end
            disp(' ')
        end

        function plot(obj,varargin)
            % PLOT method: draw the Domain object.
            %
            % domain.PLOT('labels') add nodes and edges labels in the
            % plot.
            % domain.PLOT('edgelabels') add only edges labels in the plot
            %  the plot.
            % domain.PLOT('nodelabels') add only nodes labels in the plot
            %  the plot.
            % domain.plot(I) plot only the edges specified in the index
            % vector I.



            Ne = size(obj.edges,1);

            if nargin == 1 % No input
                obj.plot(1:Ne);
                return
            end

            if ischar(varargin{1})
                varargin = [{1:Ne},varargin];
            end

            I = varargin{1};
            Ne = length(I);


            % Default parameters values
            edgelabels = false;
            nodelabels = false;
            domaincolor = [0.8 0.8 0.8];

            for n = 2 : length(varargin)
                switch varargin{n}
                    case {'edgelabels','EdgeLabels'}
                        edgelabels = true;
                    case {'nodelabels','NodeLabels'}
                        nodelabels = true;
                    case {'labels'}
                        edgelabels = true;
                        nodelabels = true;
                    case 'domaincolor'
                        if n < nv && isa(varargin{n+1},double)
                            domaincolor = varargin{n+1};
                        else
                            domaincolor = [0.2 0.2 0.8];
                        end
                end

            end


            % Plot boundary egdes
            for i = 1 : Ne
                obj.plot_edge(I(i),edgelabels,nodelabels);
                hold on
            end

            if domaincolor
                plot(obj.poly,'FaceColor',domaincolor)
            end


            axis equal
            hold off

        end
        
        function [c,obj2] = simplify(obj,h)
            obj2 = obj;
            E = obj.edges;
            E2 = {};
            c = 0;
            for i = 1 : size(E,1)
                edge = E(i,:);
                [x,y] = obj.discretize_edge(edge,h,inf);
                N = length(x);
                if N <= 3 && edge{1} == edge{2}
                    c = c + 1;
                else
                    E2 = [E2 ; edge];
                end
            end
            if c
                disp(['warning: ' num2str(c) ' small holes or inclusions have been deleted'])
                disp(['try Mesh(domain,h,h_bound) with a smaller h_bound value'])

                obj2.edges = E2;
            end

        end
        

        function [M,Ibound,obj] = geometry_matrix(obj,h)
            % GEOMTRY_MATRIX method: build the geometry matrix of the
            % domain. This matrix is needed to create a Mesh object form
            % this domain. Use: 
            %
            % [M,Ibound] = obj.MATRIX(h) compile the data contained in the Domain
            % object in a matrix M to be used by initmesh. Ibound contains
            % the indices of the i of boundary edges.

            [c,obj2] = simplify(obj,h);

            if c
                obj = obj2;
            end

            E = obj.edges;
            M = [];
            Ibound = [];
            for i = 1 : size(E,1)
                edge = E(i,:);
                side = edge{5};
                [x,y] = obj.discretize_edge(edge,h,inf);
                x = x';
                y = y';
                N = length(x);
                if N>3 
                switch side
                    case 'L'
                        sideM = [ones(1,N-1) ; zeros(1,N-1)];
                    case 'R'
                        sideM = [zeros(1,N-1) ; ones(1,N-1)];
                    case 'LR'
                        sideM = ones(2,N-1);
                end
                m = [2*ones(1,N-1) ; x(1:end-1) ; x(2:end) ;  y(1:end-1) ; y(2:end) ; sideM ; zeros(5,N-1)];
                M = [M m];
                Ibound = [Ibound ; i*ones(N-1,1)];
                end

            end

        end
    end

    methods (Access = private) % need Private !

        function obj = check_structure(obj)
            % Check the domain structure
            [Ne,dime] = size(obj.edges);
            [Nn,dim] = size(obj.nodes);
            if dime ~= 5
                error('Invalid edges format.')
            end

            if dim ~=2
                if Nn == 2
                    obj.nodes = obj.nodes';
                    Nn = dim;

                else
                    error('Invalid nodes format: first input must be a Nx2 double')
                end
            end

            for n = 1 : Ne
                e = obj.edges(n,:);
                a = min(e{1},e{2});
                b = max(e{1},e{2});
                if a<1 || b>Nn
                    error('Invalid node index in edges definiton')
                end
                if strcmp(e{3},'param')

                end

            end

        end

        function plot_edge(obj,i,edgelabels,nodelabels)
            % parameters
            Ndis = 128;


            edge = obj.edges(i,:);
            i1 = edge{1};
            i2 = edge{2};
            type = edge{3};
            side = edge{5};

            x1 = obj.nodes(i1,:);
            x2 = obj.nodes(i2,:);

            % Color choice
            switch side
                case 'L'
                    c = [0.8 0.1 0.1];

                case 'R'
                    c = [0.1 0.8 0.1];

                case 'LR'
                    c = [0.1 0.1 0.8];
                otherwise
                    error(['Wrong last element of the edge cell ' num2str(i) ]);
            end

            % Edge orientation
            a = angle((x2(1) - x1(1)) +1i*(x2(2) - x1(2)));
            ha = 'center';
            va = 'bottom';
            if abs(a)<pi/4
                ha = 'center';
                va = 'top';
            end
            if a>= pi/4 && a<3*pi/4
                ha = 'left';
                va = 'middle';
            end
            if a<= -pi/4 && a>-3*pi/4
                ha = 'right';
                va = 'middle';
            end

            switch type


                case 'straight'

                    plot([x1(1) x2(1)],[x1(2) x2(2)],'Color',c,'linewidth',2);
                    hold on
                    plot([x1(1) x2(1)],[x1(2) x2(2)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');

                    if edgelabels
                        text((x1(1)+x2(1))/2,(x1(2)+x2(2))/2,[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                    end
                    if nodelabels
                        text(x1(1),x1(2),[ ' ' num2str(i1) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');
                    end


                case 'param'

                    [x,y] = discretize_edge(obj,edge,0,Ndis);
                    plot(x,y,'Color',c,'linewidth',2)
                    hold on

                    plot([x(1) x(end)],[y(1) y(end)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if edgelabels
                        text(x(Ndis/2),y(Ndis/2),[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                    end
                    if nodelabels
                        text(x1(1),x1(2),[ ' ' num2str(i1) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');
                    end
                case 'data'
                    data = edge{4};
                    Ndata = size(data,1);
                    x = data(:,1);
                    y = data(:,2);
                    plot(x,y,'Color',c,'linewidth',2)
                    hold on

                    plot([x(1) x(end)],[y(1) y(end)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if edgelabels
                        text(x(floor(Ndata/2)),y(floor(Ndata/2)),[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                    end
                    if nodelabels
                        text(x1(1),x1(2),[ ' ' num2str(i1) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');
                    end

            end
        end

        function [x,y] = discretize_edge(obj,edge,h,Ndis)
            type = edge{3};
            param = edge{4};
            
            switch type
                case 'straight'
                    x1 = obj.nodes(edge{1},:);
                    x2 = obj.nodes(edge{2},:);
                    l = sqrt((x2(1)-x1(1)).^2 + (x2(2)-x1(2)).^2);
                    N = max(floor(l/h)+1,2);
                    N = min(N,Ndis);
                    x = linspace(x1(1),x2(1),N)';
                    y = linspace(x1(2),x2(2),N)';
                case 'param'
                    interval = param.interval;
                    tmin = interval(1);
                    tmax = interval(2);
                    if Ndis == inf
                        Ndis = 256;
                    end
                    t = linspace(tmin,tmax,Ndis)';
                    [teq,L,x,y] = iterative_equalizer(t,param);
                    if h
                        l = sum(L);
                        N = 1 + ceil(l/h);
                        teq = interp1(t,teq,linspace(tmin,tmax,N)');
                        [teq,L,x,y] = iterative_equalizer(teq,param);
                    end
                case 'data'
                    x = param(:,1);
                    y = param(:,2);
                    L = sqrt(diff(x).^2+diff(y).^2);
                    if min(L)<h
                        p = under_sample2(param,h);
                        x = p(:,1);
                        y = p(:,2);
                    end
            end

        end

        function P = build_poly(obj)
            warning('off');
            Ndis = 256;
            Ne = size(obj.edges,1);
            Px = {};
            Py = {};
            pol = [];
            for n = 1 : Ne
                edge = obj.edges(n,:);
                type = edge{3};
                side = edge{5};
                if ~strcmp(side,'LR')
                    i = edge{1};
                    j = edge{2};
                    pi = obj.nodes(i,:);
                    pj = obj.nodes(j,:);
                    if isempty(pol)
                        istart = i;
                        pol = pi;
                    end

                    switch type
                        case 'param'
                            [x,y] = obj.discretize_edge(edge,0,Ndis);
                            pol = [pol ; [x(2:end-1) y(2:end-1)]];
                        case 'data'
                            data = edge{4};
                            x = data(:,1);
                            y = data(:,2);
                            pol = [pol ; [x(2:end-1) y(2:end-1)]];
                    end



                    if j~=istart
                        pol = [pol ; pj];
                    else
                        Px = [Px pol(:,1)'];
                        Py = [Py pol(:,2)'];
                        pol = [];
                    end
                end

            end
            P = polyshape(Px,Py);
            warning('on');
        end

    end
end



% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ef = format_edge(e)
ef = cell(1,5);
ne = length(e);
ef{1} = e{1};
ef{2} = e{2};
ef{5} = 'L';
last = e{ne};
if ischar(last)
    ne = ne-1;
    switch last
        case {'l','left','L'}

        case {'r','right','R'}
            ef{5} = 'R';
        case {'rl','lr','LR','RL','leftright','rightleft','both'}
            ef{5} = 'LR';
        otherwise
            error('Wrong last element of the edge cell.')
    end
end

if ne == 2
    ef{3} = 'straight';
    return
end

w = class(e{3});
switch w
    case 'cell'
        e3cell = e{3};
        param.x = e3cell{1};
        param.y = e3cell{2};
        param.interval = e{4};
        ef{3} = 'param';
        ef{4} = param;
    case 'char'
        param.x = e{3};
        param.y = e{4};
        param.interval = e{5};
        ef{3} = 'param';
        ef{4} = param;

    case 'function_handle'
        if ne == 4
            param.xy = e{3};
            param.interval = e{4};
        elseif ne == 5
            param.x = e{3};
            param.y = e{4};
            param.interval = e{5};
        else
            error('Wrong 3rd element of the edge cell.')
        end
        ef{3} = 'param';
        ef{4} = param;
    case 'double'
        ef{3} = 'data';
        ef{4} = e{3};

    otherwise
        error('Wrong 3rd element of the edge cell.')
end
end

function [x,y] = eval_edge(param,t)
if isfield(param,'xy')
    p = feval(param.xy,t);
    x = p(:,1);
    y = p(:,2);
    return
end

expx = param.x;
expy = param.y;
wx = whos('expx');
if strcmp(wx.class,'char')
    x = eval(expx);
elseif strcmp(wx.class,'function_handle')
    x = feval(expx,t);
end
wy = whos('expy');
if strcmp(wy.class,'char')
    y = eval(expy);
elseif strcmp(wy.class,'function_handle')
    y = feval(expy,t);
end
end

function [teq,L,x,y] = iterative_equalizer(t,param)
[x,y] = eval_edge(param,t);
L = sqrt(diff(x).^2+diff(y).^2);
teq = t;
n = 0;
while max(L)/min(L) > 1 + 1e-3
    teq = equalizer(teq,x,y);
    [x,y] = eval_edge(param,teq);
    L = sqrt(diff(x).^2+diff(y).^2);
    n = n+1;
    if n == 5
        disp('Warning: max iteration reached')
        break
    end
end
end

function teq = equalizer(t,x,y)
N = length(t);
l_list = sqrt(diff(x).^2+diff(y).^2);
L = sum(l_list);
dx0 = L/(N-1);
teq = t;
i = 1;
n = 2;
dx = dx0;
l = l_list(1);
tau = t(1);
while n < N
    if dx<l
        teq(n) = tau + dx/l*(t(i+1)-tau);
        tau = teq(n);
        l = l - dx;
        dx = dx0;
        n = n + 1;
    else
        dx = dx - l;
        i = i + 1;
        l = l_list(i);
        tau = t(i);
    end
end
end

function ph = under_sample(p,h)
if size(p,1) == 0
    ph = [];
    return
end
ph = p(1,:);
P = p(2:end,:);
while size(P,1)>0
    p1 = ph(end,:);
    p2 = P(1,:);
    d = norm((p2-p1));
    if d >= h
        ph = [ph ; p2];
    end
    P = P(2:end,:);
end

ph(1,:) = (ph(1,:) + ph(end,:))/2;
ph = ph(1:end-1,:);
end

function ph = under_sample2(p,h)
N = size(p,1);
ph = p(1,:);
x = p(1,:);
d = 0;
n = 1;
while n<N
    while d<h && n<N
        n = n + 1;
        z = p(n,:);
        d = norm(x - z);
    end
    if n == N
        if d<h/2
        ph = [ph(1:end-1,:); p(n,:)];
        else
            ph = [ph ; p(n,:)];
        end

    else
        y = p(n-1,:);
        v = (z - y)/norm(z - y);
        px = y + ((x - y)*v')*v;
        alpha = sqrt(h^2 - norm(px - x)^2);
        x_new = px + alpha*v;
        ph = [ph ; x_new];
        x = x_new;
        d = 0;
    end
end


end