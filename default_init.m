function [contract,parameter,method,grid] = default_init()
%Constructor that initiates default problem setup
%   contract: contract
%   parameter: parameter
%   method: method
%   grid: grid

%contract
contract.payoff = 'EUcall';
contract.K = 100;
contract.T = 1;

%parameters
parameter.sig = 0.3;
parameter.r = 0.01;

%method
method.space.type = 'FD';
method.space.order = 2;
method.time.type = 'BDF';
method.time.order = 1;
method.time.M = 10;

%grid
grid.type = 'uniform';
grid.dim = 1;
grid.N = 101;
grid.xmax = 4*contract.K;
end

