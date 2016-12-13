function [grid] = make_grid(grid)
%Creates a computational grid.
%Inputs: grid
%   grid.dim: 1, 2, 3;
%   grid.type: 'uniform';
%   grid.N: number of points;
%   grid.xmax: far field boundary.
%Outputs:
%   grid.x: node coordinates;
%   grid.ind: indices of all nodes;
%   grid.indin: indices of nodes inside the domain;
%   grid.indcf: indices of nodes at close field boundary;
%   grid.indff: indices of nodes at far field boundary.

%% Grid

switch grid.dim
    case 1
        switch grid.type
            case 'uniform'
                grid.x = transpose(linspace(0,grid.xmax,grid.N));
                grid.ind = 1:grid.N;
                grid.indin = 2:(grid.N-1);
                grid.indcf = 1;
                grid.indff = grid.N;
                grid.Ntot = grid.N;
                grid.h = grid.x(2)-grid.x(1);
                return;
        end
        
        
        
    case 2
        switch grid.type
            case 'uniform'
                x = linspace(0,grid.xmax,grid.N);
                y = linspace(0,grid.xmax,grid.N);
                
                [X,Y] = meshgrid(x,y);
                xvec = X(:);
                yvec = Y(:);
                
                grid.ind = 1:numel(xvec);
                
%                 indle = [1:grid.N];
%                 inddo = [1:grid.N:grid.N^2-grid.N+1];
                indup = grid.N:grid.N:grid.N^2;
                indri = grid.N^2-grid.N+1:1:grid.N^2;
                
                grid.indcf = 1;
                grid.indff = [indup(1:end-1),indri];
                
                grid.indin = grid.ind; grid.indin([grid.indff,grid.indcf]) = [];
                
                grid.Ntot = numel(xvec);
                grid.x = [xvec,yvec];
                grid.h = [grid.x(2,1) - grid.x(1,1), grid.x(2,2) - grid.x(1,2)];
                return;
        end
        
        
        
    case 3
        disp('Not implemented.');
        return;
end
end
