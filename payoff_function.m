function [u] = payoff_function(payoff,x)
%Creates a computational grid.
%Inputs:
%   payoff: 'EUcall';
%   x: node coordinates returned by make_grid;
%Outputs:
%   u: values of payoff at node coordinates x;

switch payoff
    case 'EUcall'
        u = max(x-K,0);
end
end
