function [u] = payoff_function(c,s)
%Sets the terminal condition on the grid.
%Inputs: c, s
%   c.payoff: 'EUcall'
%   c.K: strike
%   c.T: maturity
%   x: node coordinates returned by make_grid;
%Outputs:
%   u: values of payoff at node coordinates x;

switch c.payoff
    
    case {'EUcall', 'AMcall'}
        u = max(s.x - c.K, 0);
        return;
        
    case {'EUCallbasket', 'AMcallbasket'}
        u = max((1/s.dim)*sum(s.x,2) - c.K, 0);
        return;
end
end
