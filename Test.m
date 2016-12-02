clc
clear
close all

%% Test EUcall

%#% Setup
%contract
c.payoff = 'EUcall';
c.K = 100;
c.T = 1;

%parameters
p.sig = 0.3;
p.r = 0.01;

%grid
g.type = 'uniform';
g.dim = 1;
g.N = 100;
g.smax = 4*c.K;

%method
m.type = 'FD';
m.order = 2;


%#% Evaluation
s = make_grid(g);
s.u = payoff_function(c,s);
s.W = differentiation_matrix(m,p,s);
