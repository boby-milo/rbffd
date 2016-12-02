function [s] = make_grid(g)
%Creates a computational grid.
%Inputs: g
%   g.dim: 1, 2, 3;
%   g.type: 'uniform';
%   g.N: number of points;
%   g.smax: far field boundary.
%Outputs:
%   s.x: node coordinates;
%   s.ind: indices of all nodes;
%   s.indin: indices of nodes inside the domain;
%   s.indcf: indices of nodes at close field boundary;
%   s.indff: indices of nodes at far field boundary.

%% Grid
switch g.dim
    
    case 1
        switch g.type
            case 'uniform'
                s.x = transpose(linspace(0,g.smax,g.N));
                s.ind = 1:g.N;
                s.indin = 2:(g.N-1);
                s.indcf = 1;
                s.indff = g.N;
                s.type = g.type;
                s.N = g.N;
                s.h = s.x(2)-s.x(1);
                return;
        end
        
        
        
    case 2
        switch g.type
            case 'uniform'
                x = linspace(0,g.smax,g.N);
                y = linspace(0,g.smax,g.N);
                
                [X,Y] = meshgrid(x,y);
                xvec = X(:);
                yvec = Y(:);
                
                s.ind = 1:numel(xvec);
                
%                 indle = [1:g.N];
%                 inddo = [1:g.N:g.N^2-g.N+1];
                indup = [g.N:g.N:g.N^2];
                indri = [g.N^2-g.N+1:1:g.N^2];
                
                s.indcf = 1;
                s.indff = [indup(1:end-1),indri];
                
                s.indin = s.ind; s.indin([s.indff,s.indcf]) = [];
                
                s.N = numel(xvec);
                s.x = [xvec,yvec];
                return;
        end
end
end
