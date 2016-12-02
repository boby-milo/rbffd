function [x,ind,indin,indcf,indff] = make_grid(dim,type,N,smax)
%Creates a computational grid.
%Inputs:
%   dim: 1, 2, 3;
%   type: 'uniform';
%   N: number of points;
%   smax: far field boundary.
%Outputs:
%   x: node coordinates;
%   ind: indices of all nodes;
%   indin: indices of nodes inside the domain;
%   indcf: indices of nodes at close field boundary;
%   indff: indices of nodes at far field boundary.

%% Grid
switch dim
    case 1
        switch type
            case 'uniform'
                x = transpose(linspace(0,smax,N));
                ind = 1:N;
                indin = 2:N-1;
                indcf = 1;
                indff = N;
        end
end
end

