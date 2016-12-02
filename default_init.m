function [c,p,m,s] = default_init()
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
m.time.M = 10;

%grid
s.type = 'uniform';
s.dim = 1;
s.N = 100;
s.xmax = 4*c.K;
end

