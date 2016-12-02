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
end
end
