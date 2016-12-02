clc
clear
close all

dbstop if error

%% Test EUcall

%#% Setup
[c,p,m,g] = default_init();

%#% Evaluation
s = make_grid(g);
s.u0 = payoff_function(c,s);
s.W = differentiation_matrix(m,p,s);
[s.u,m] = time_integrate(m,c,p,s);
