function [s] = make_grid(s)
%Creates a computational grid.
%Inputs: s
%   s.dim: 1, 2, 3;
%   s.type: 'uniform';
%   s.N: number of points;
%   s.smax: far field boundary.
%Outputs:
%   s.x: node coordinates;
%   s.ind: indices of all nodes;
%   s.indin: indices of nodes inside the domain;
%   s.indcf: indices of nodes at close field boundary;
%   s.indff: indices of nodes at far field boundary.

%% Grid

switch s.dim
    case 1
        switch s.type
            case 'uniform'
                s.x = transpose(linspace(0,s.smax,s.N));
                s.ind = 1:s.N;
                s.indin = 2:(s.N-1);
                s.indcf = 1;
                s.indff = s.N;
                s.N = s.N;
                s.h = s.x(2)-s.x(1);
                return;
        end
        
        
        
    case 2
        switch s.type
            case 'uniform'
                x = linspace(0,s.smax,s.N);
                y = linspace(0,s.smax,s.N);
                
                [X,Y] = meshgrid(x,y);
                xvec = X(:);
                yvec = Y(:);
                
                s.ind = 1:numel(xvec);
                
%                 indle = [1:s.N];
%                 inddo = [1:s.N:s.N^2-s.N+1];
                indup = [s.N:s.N:s.N^2];
                indri = [s.N^2-s.N+1:1:s.N^2];
                
                s.indcf = 1;
                s.indff = [indup(1:end-1),indri];
                
                s.indin = s.ind; s.indin([s.indff,s.indcf]) = [];
                
                s.Ntot = numel(xvec);
                s.x = [xvec,yvec];
                s.h = [s.x(2,1) - s.x(1,1), s.x(2,2) - s.x(1,2)];
                return;
        end
end
end
