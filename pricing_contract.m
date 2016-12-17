function [u] = pricing_contract(contract,grid, smoothing)
%Sets the terminal condition on the grid. This works for n-dimensional
%basket options.
%Inputs: contract, grid
%   contract.payoff: 'EUcall'
%   contract.K: strike
%   contract.T: maturity
%   x: node coordinates returned by pricing_grid;
%Outputs:
%   u: values of payoff at node coordinates x;

if nargin <= 2
    smoothing = 0;
end

switch contract.payoff
    case {'EUcall', 'AMcall'}
        f = @(x) max(x - contract.K,0);

    case {'EUput', 'AMput'}
        f = @(x) max(contract.K - x,0);

    case {'EUcallBasket', 'AMcallBasket'}
        u = max((1/grid.dim)*sum(grid.x,2) - contract.K, 0);
        return;

    case {'EUputBasket', 'AMputBasket'}
        u = max(contract.K - (1/grid.dim)*sum(grid.x,2), 0);
        return;
end


if smoothing
    if grid.dim == 1
        a = 1/10;
        indreg = find((grid.x >= (1-a)*contract.K) & (grid.x <= (1+a)*contract.K));
        xind = grid.x(indreg);
        uind = pricing_grid_smooth4(xind,f);
        u = f(grid.x);
        u(indreg)=uind;
        return;
    else
        disp('Higher-dimensional smoothing is not implemented.')
        return;
    end

else
    u = f(grid.x);
    return;
end

end
