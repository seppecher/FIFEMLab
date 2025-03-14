classdef Domain

    properties (SetAccess = private)
        nodes
        edges
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
                obj.nodes  = [0 0 ; 1 0 ; 1 1 ; 0 1];
                obj.edges = {{1,2};{2,3};{3,4};{4,1}};
                return
            end

            arg1 = varargin{1};

            if ischar(arg1)
                if nargin == 1
                    obj = Domain(arg1,1);
                    return
                end
                r = varargin{2};
                switch arg1
                    case {'square','Square'}
                        obj.nodes  = [0 0 ; r 0 ; r r ; 0 r];
                        obj.edges = {{1,2};{2,3};{3,4};{4,1}};
                        return
                    case {'disc','Disc'}
                            edgefun = @(t) r*[cos(t) sin(t)];
                            obj.nodes  = [r 0];
                            obj.edges = {{1,1,edgefun,[0,2*pi]}};
                        return
                    otherwise
                        error('Unvalid first argument.')
                end

            end

            if nargin == 1
                nodes = arg1;
                if size(nodes,2) ~= 2
                    nodes = nodes';
                    if size(nodes,2) ~=2
                        error('Wrong format for first input: require a Nx2 nodes matrix.');
                    end
                end
                Nn = size(nodes,1);
                edges = cell(Nn,1);
                for n = 1 : Nn-1
                    edges{n} = {n,n+1};
                end
                edges{Nn} = {Nn,1};
                obj.nodes = nodes;
                obj.edges = edges;
                return
          

            end




            % Basic constructor
            obj.nodes = arg1;
            obj.edges = varargin{2};
            obj = obj.check_structure;

        end

        function disp(obj)
            % DOMAIN.DISP: display informations contained in the Domain
            % object. 
            disp('Domain object')

            Ne = length(obj.edges);
            if Ne <= 50 % Full description
                


                disp(['Nodes (' num2str(size(obj.nodes,1)) '):']);
                disp(obj.nodes)
                disp(['Edges (' num2str(size(obj.edges,1)) '):']);
                for i = 1 : length(obj.edges)
                    edge = obj.edges{i};
                    last = edge{end};
                    if ischar(last)
                        edge = edge(1:end-1);
                    end
                    a = num2str(edge{1});
                    b = num2str(edge{2});
                    type = '';
                    if ischar(last)
                       
                        switch last
                            case {'r','right','R'}
                                type = ', type = hole';
                            case {'rl','lr','LR','RL','leftright','rightleft','both'}
                                type = ', type = internal';
                        end
                    end
                    switch length(edge)
                        case {2}
                            disp(['edge ' num2str(i) ': from ' a ' to ' b ' straight' type])
                        case {3}
                            disp(['edge ' num2str(i) ': from ' a ' to ' b ' curve data' type])
                        case {4,5,6}
                            disp(['edge ' num2str(i) ': from ' a ' to ' b ' parametric' type])
                    end

                end
            else % Synthetic description
                disp(['Nodes number = ' num2str(size(obj.nodes,1))]);
                disp(['Edges number = ' num2str(size(obj.edges,1))]);
            end
            disp(' ')
        end

        function plot(obj,varargin)
            % PLOT method: draw the Domain object. 
            % 
            % domain.PLOT('labels') add nodes and edges labels in the
            % plot.
            % domain.PLOT('edgeslabels') add only edges labels in the plot
            %  the plot.
            % domain.PLOT('nodeslabels') add only nodes labels in the plot
            %  the plot.
            % domain.plot(I) plot only the edges specified in the index
            % vector I.

            Ne = length(obj.edges);

            if nargin == 1 % No input
                obj.plot(1:Ne);
                return
            end
            if ischar(varargin{1}) 
                obj.plot(1:Ne,varargin);
                return
            end

            I = varargin{1};
            Ne = length(I);

            for i = 1 : Ne
                obj.plot_edge(I(i),varargin{2:end});
                hold on
            end
            axis equal
            hold off

        end
 
    end
    
    methods (Access = private)

        function obj = check_structure(obj)

            % Check the domain structure
            Ne = length(obj.edges);
            [Nn,dim] = size(obj.nodes);
            if dim ~=2
                if Nn == 2
                    obj.nodes = obj.nodes';
                    Nn = dim;

                else
                    error('Invalid node format: first input must be a Nx2 double')
                end
            end

            for n = 1 : Ne
                e = obj.edges{n};
                a = min(e{1},e{2});
                b = max(e{1},e{2});
                if a<1 || b>Nn
                    error('Invalid node index in edges definiton')
                end
            end

        end



        function plot_edge(obj,i,label)
            edge = obj.edges{i};
            le = length(edge);
            last = edge{end};

            % Color choice
            if ischar(last)
                switch last
                    case {'l','left','L'}
                            c = [0.8 0.1 0.1];
                        case {'r','right','R'}
                            c = [0.1 0.8 0.1];
                        case {'rl','lr','LR','RL','leftright','rightleft','both'}
                            c = [0.1 0.1 0.8];
                    otherwise
                        error(['Wrong last element of the edge cell ' num2str(i) ])
                end
            else 
                c = [0.8 0.1 0.1];
            end

            x1 = obj.nodes(edge{1},:);
            x2 = obj.nodes(edge{2},:);
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


            switch le
                
                case {2,3}
                    
                    plot([x1(1) x2(1)],[x1(2) x2(2)],'Color',c,'linewidth',2);
                     hold on
                    plot([x1(1) x2(1)],[x1(2) x2(2)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if nargin == 3
                        
                        switch label{1} 

                            case 'edgelabels'
                                text((x1(1)+x2(1))/2,(x1(2)+x2(2))/2,[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                            case 'nodelabels'
                                text(x1(1),x1(2),[ ' ' num2str(edge{1}) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');
                            otherwise
                                 text((x1(1)+x2(1))/2,(x1(2)+x2(2))/2,[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                                 text(x1(1),x1(2),[ ' ' num2str(edge{1}) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');    
                        end
                    end
                    
                case {4,5,6}
                    Ndis = 256;
                   [x,y] = discretize_edge(obj,edge,0,Ndis);
                    plot(x,y,'Color',c,'linewidth',2)
                    hold on
                    
                    plot([x(1) x(end)],[y(1) y(end)],'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k');
                    if nargin == 3
                        switch label{1} 

                            case 'edgelabels'
                                text(x(Ndis/2),y(Ndis/2),[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                            case 'nodelabels'
                                text(x1(1),x1(2),[ ' ' num2str(edge{1}) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va,'fontweight','bold');
                            otherwise
                                text(x(Ndis/2),y(Ndis/2),[ ' ' num2str(i) ' ' ],'Color',c,'FontSize',12,'HorizontalAlignment',ha,'VerticalAlignment',va);
                                text(x1(1),x1(2),[ ' ' num2str(edge{1}) ' ' ],'Color','k','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle','fontweight','bold');    
                        end
                        
                    end
            end

        end


        function [x,y] = discretize_edge(obj,edge,dx,Ndis)
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
            N = Ndis;
            t = linspace(tmin,tmax,N)';
            [x,y] = eval_edge(edgedef,t);
            tmin = tdom(1);
            tmax = tdom(2);
            
            h = (1e-10)*(tmax-tmin);
            th = t+h;
            [xh,yh] = eval_edge(edgedef,th);
            
            xp = (xh - x)/h;
            yp = (yh - y)/h;
            
            v = sqrt(xp.^2+yp.^2);
            v = (v(2:end) + v(1:end-1))/2;
            l = sum(v)/(N-1)*(tmax-tmin);
            if dx == 0
                N = Ndis;
            else
                N = max(floor(l/dx)+2,3);
            end
            t = linspace(tmin,tmax,N)';
            [x,y] = eval_edge(edgedef,t);
            h = (1e-10)*(tmax-tmin);
            th = t+h;
            [xh,yh] = eval_edge(edgedef,th);
            
            xp = (xh - x)/h;
            yp = (yh - y)/h;
            
            v = sqrt(xp.^2+yp.^2);
            v = (v(2:end) + v(1:end-1)/2);
            cv =[0 ; cumsum(1./v)];
            cv = cv/max(cv);
            t = (cv*(tmax-tmin)+tmin);
            
            [x,y] = eval_edge(edgedef,t);
            
        end

    end
end



% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = eval_edge(edgedef,t)
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