function [c,p,m,g] = default_init()
%Constructor that initiates default problem setup
%   c: contract
%   p: parameters
%   m: method
%   g: grid

%contract
c.payoff = 'EUcall';
c.K = 100;
c.T = 1;

%parameters
p.sig = 0.3;
p.r = 0.01;

%method
m.space.type = 'FD';
m.space.order = 2;
m.time.type = 'BDF';
m.time.order = 1;
m.time.M = 100;

%grid
g.type = 'uniform';
g.dim = 1;
g.N = 1000;
g.smax = 4*c.K;
end

