classdef Domain

    properties (SetAccess = private)
        nodes
        edges
    end
    
    methods (Access = public)
        function obj = Domain(input1,input2)
            % Class DOMAIN constructor. This class allows to simply
            % describe various 2D domains. A Domain object is needed to create
            % later a Mesh object. 
            %
            % PROPERTIES: 
            % nodes
            % edges
            %
            % PUBLIC METHODS: 
            % disp
            % plot
            %
            % BASIC USE OF THE CONSTRUCTOR:
            % d = Domain(nodes,edges) where nodes is a Nx2 matrix of N
            % nodes positions and edges is a Nx2 matrix of indices of nodes
            % defining the edges in trigonometric order.
            %
            % Run script sample_domain.m and read fifemtoolbox
            % documentation for a complete descrition

           
            if nargin == 0
                obj = Domain('square');
                return
            end
            w = whos('input1');
            switch w.class
                case 'char'
                    switch input1
                        case {'square','Square'}
                            if nargin == 2
                                r = input2;
                            else
                                r = 1;
                            end
                            input1 = [-r -r;r -r;r r;-r r];
                            input2 = {{1,2};{2,3};{3,4};{4,1}};
                            obj = Domain(input1,input2);
                            
                        case {'disc','Disc'}
                            if nargin == 2
                                r = input2;
                            else
                                r = 1;
                            end
                            input1 = [r 0;];
                            edgefun = @(t) r*[cos(t) sin(t)];
                            input2 = {{1,1,edgefun,[0,2*pi]}};
                            obj = Domain(input1,input2);
                    end
                    
                case 'double'
                    switch nargin
                        case 1
                            obj.nodes = input1;
                            N =  size(input1,1);
                            edge = cell(N,1);
                            for i = 1 : N-1
                                edge{i} = {i,i+1};
                            end
                            edge{N} = {N,1};
                            obj.edges = edge;
                            
                        case 2
                            obj.nodes = input1;

                            w2 = whos('input2');
                            switch w2.class
                                case 'double'
                                    n = size(input2,1);
                                    if n == 2
                                        input2 = input2';
                                        n = p;
                                    end

                                    C = cell(1,n);
                                    for k = 1 : n
                                        C{k} = {input2(k,1),input2(k,2)};
                                    end
                                     obj.edges = C;
                                case 'cell'
                                    n = size(input2,1);
                                    if n>1
                                        input2 = input2';
                                    end
                                    obj.edges = input2;
                                otherwise
                                    error('Wrong type for input 2.');
                            end
                    end
                    
                otherwise
                    error('First input class must be double or char')
            end
            
            % Validity test
            Ne = length(obj.edges);
            Nn = size(obj.nodes,1);
            for n = 1 : Ne
                e = obj.edges{n};
                a = min(e{1},e{2});
                b = max(e{1},e{2});
                if a<1 || b>Nn
                    error('Invalid node index')
                end
            end
            
            
        end
        function disp(obj)
            % DOMAIN.DISP:
            %       obj.disp display informations contained in the Domain
            %       object.
            disp('Domain object')
            disp('Nodes:')
            disp(obj.nodes)
            for i = 1 : length(obj.edges)
                edge = obj.edges{i};
                a = num2str(edge{1});
                b = num2str(edge{2});
                switch length(edge)
                    case {2,3}
                        disp(['edge ' num2str(i) ': from ' a ' to ' b])
                    case {4,5,6}
                        disp(['edge ' num2str(i) ': from ' a ' to ' b ' parametric'])
                end
                
            end
            disp(' ')
        end
        function plot(obj,arg2,arg3)
             % DOMAIN.PLOT:
            %       obj.PLOT plot the Domain object. 
            %       obj.PLOT('labels') plot the Domain object with edges
            %       labels.
            switch nargin
                case 1
                    obj.plot(1:length(obj.edges));
                case 2
                    if strcmp(arg2,'labels')
                        obj.plot(1:length(obj.edges),'labels');
                        return
                    end
                    I = arg2;
                    Ne = length(I);
                    for i = 1 : Ne
                        obj.plotedge(obj.edges{I(i)});
                        hold on
                    end
                    axis equal
                    hold off
                case 3
                    if strcmp(arg3,'labels')
                        I = arg2;
                        Ne = length(I);
                        for i = 1 : Ne
                            obj.plotedge(obj.edges{I(i)},num2str(I(i)));
                            hold on
                        end
                        axis equal
                        hold off
                        return
                    else
                        error('Wrong argument')
                    end
                    axis equal
                    hold off       
                    
            end
            
        end
        function [M,Ibound]=matrix(obj,dx)
            % DOMAIN.MATRIX:
            %       [M,Ibound]=obj.MATRIX compile the data contained in the Domain
            %       object in a matrix M to be used by initmesh. Ibound is
            %       the number of the initial edge.
            E = obj.edges;
            M = [];
            Ibound = [];
            for i = 1 : length(E)
                edge = E{i};
                [x,y] = obj.discretize_edge(edge,dx);
                x = x';
                y = y';
                N = length(x);
                
                last = edge{end};
                if ischar(last)
                    switch last
                        case {'l','left','L'}
                            side = [ones(1,N-1) ; zeros(1,N-1)];
                        case {'r','right','R'}
                            side = [zeros(1,N-1) ; ones(1,N-1)];
                        case {'rl','lr','LR','RL','leftright','rightleft'}
                            side = ones(2,N-1);
                        otherwise
                            error('Wrong last argument')
                    end
                    m = [2*ones(1,N-1) ; x(1:end-1) ; x(2:end) ;  y(1:end-1) ; y(2:end) ; side ; zeros(5,N-1)];
                else
                    m = [2*ones(1,N-1) ; x(1:end-1) ; x(2:end) ;  y(1:end-1) ; y(2:end) ; ones(1,N-1) ; zeros(6,N-1)];
                end
                M = [M m];
                Ibound = [Ibound ; i*ones(N-1,1)];
            end
        end
        function b = eq(obj,domain)
            if size(obj.nodes,1) == size(domain.nodes,1)
                s = sum(sum(abs(obj.nodes-domain.nodes)));
            else
                s = 1;
            end
            b = ~(s>0);
            e = obj.edges;
            if b
                for i = 1 : length(e)
                    a = (e{i}{1} ~= domain.edges{i}{1}) || (e{i}{2} ~= domain.edges{i}{2});
                    if a
                        b=0;
                        break
                    end
                end
            end
        end
    end
    
    methods (Access = private)
        function [x,y] = evaledge(obj,edgedef,t)
            w = whos('edgedef');
            if strcmp(w.class,'cell')
                expx = edgedef{1};
                expy = edgedef{2};
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
            elseif strcmp(w.class,'function_handle')
                p = feval(edgedef,t);
                x = p(:,1);
                y = p(:,2);
            end
            
        end
        function plotedge(obj,edge,label)
            le = length(edge);
            last = edge{end};
            if ischar(last)
                switch last
                    case {'l','left','L'}
                            c = [0.8 0 0];
                        case {'r','right','R'}
                            c = [0 0.8 0];
                        case {'rl','lr','LR','RL','leftright','rightleft'}
                            c = [0 0 0.8];
                    otherwise
                        error('Wrong last argument')
                end
            else 
                c = [0.8 0 0];
            end
                    
            switch le
                
                case {2,3}
                    x1 = obj.nodes(edge{1},:);
                    x2 = obj.nodes(edge{2},:);
                    plot([x1(1) x2(1)],[x1(2) x2(2)],'Color',c,'linewidth',2);
                    hold on
                    plot([x1(1) x2(1)],[x1(2) x2(2)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if nargin == 3
                            text((x1(1)+x2(1))/2,(x1(2)+x2(2))/2,[' ' label],'Color',c,'FontSize',12,'VerticalAlignment','bottom');
                    end
                case {4,5,6}
                   [x,y] = discretize_edge(obj,edge,0);
                    plot(x,y,'Color',c,'linewidth',2)
                    hold on
                    
                    plot([x(1) x(end)],[y(1) y(end)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if nargin == 3
                        text(x(250),y(250),[' ' label],'Color',c,'FontSize',12,'VerticalAlignment','bottom');
                    end
            end
            
            
            
        end
        function [x,y] = discretize_edge(obj,edge,dx)
            if ischar(edge{end})
                edge = edge(1:end-1);
            end
            le = length(edge);
            switch le
                case 2
                    x1 = obj.nodes(edge{1},:);
                    x2 = obj.nodes(edge{2},:);
                    l = sqrt((x2(1)-x1(1)).^2 + (x2(2)-x1(2)).^2);
                    N = max(floor(l/dx)+1,2);
                    x = linspace(x1(1),x2(1),N)';
                    y = linspace(x1(2),x2(2),N)';
                    return
                    
                case 4
                    edgedef = edge{3};
                    tdom = edge{4};
                case 5
                    edgedef = {edge{3},edge{4}};
                    tdom = edge{5};
                otherwise
                    disp(edge)
                    error('Incorect edge length (cell of maximum length 7)')
            end
            tmin = tdom(1);
            tmax = tdom(2);
            N = 500;
            t = linspace(tmin,tmax,N)';
            [x,y] = obj.evaledge(edgedef,t);
            tmin = tdom(1);
            tmax = tdom(2);
            
            h = (1e-10)*(tmax-tmin);
            th = t+h;
            [xh,yh] = obj.evaledge(edgedef,th);
            
            xp = (xh - x)/h;
            yp = (yh - y)/h;
            
            v = sqrt(xp.^2+yp.^2);
            v = (v(2:end) + v(1:end-1))/2;
            l = sum(v)/(N-1)*(tmax-tmin);
            if dx == 0
                N = 500;
            else
                N = max(floor(l/dx)+2,3);
            end
            t = linspace(tmin,tmax,N)';
            [x,y] = obj.evaledge(edgedef,t);
            h = (1e-10)*(tmax-tmin);
            th = t+h;
            [xh,yh] = obj.evaledge(edgedef,th);
            
            xp = (xh - x)/h;
            yp = (yh - y)/h;
            
            v = sqrt(xp.^2+yp.^2);
            v = (v(2:end) + v(1:end-1)/2);
            cv =[0 ; cumsum(1./v)];
            cv = cv/max(cv);
            t = (cv*(tmax-tmin)+tmin);
            
            [x,y] = obj.evaledge(edgedef,t);
            
        end
    end
end


