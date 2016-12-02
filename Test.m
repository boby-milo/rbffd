clc
clear
close all

dbstop if error

%% Test EUcall

%#% Setup
[c,p,m,s] = default_init();

%#% Evaluation
s = make_grid(s);
s.u0 = payoff_function(c,s);
s.W = differentiation_matrix(m,p,s);
[s.u,m] = time_integrate(m,c,p,s);

figure()
plot(s.x,s.u0,s.x,s.u)

disp('EUcall test passed!')

%% Test EUbasket

%#% Setup
[c,p,m,s] = default_init();
c.payoff = 'EUbasket';

p.sig = [0.3, 0.3];
p.rho = 0.5;

s.dim = 2;
s.smax = 8*c.K;

%#% Evaluation
s = make_grid(s);
s.u0 = payoff_function(c,s);
s.W = differentiation_matrix(m,p,s);
% [s.u,m] = time_integrate(m,c,p,s);
disp('EUbasket test passed!')

figure()
clf
plot(s.x(s.indff,1),s.x(s.indff,2),'b^','MarkerFaceColor','auto')
hold on
plot(s.x(s.indcf,1),s.x(s.indcf,2),'rsq','MarkerFaceColor','auto')
plot(s.x(s.indin,1),s.x(s.indin,2),'ko','MarkerFaceColor','auto')
plot(s.x(:,1),s.x(:,2),'k.')
axis equal
axis tight
hold off

figure()
tri = delaunay(s.x(:,1),s.x(:,2));
trisurf(tri, s.x(:,1),s.x(:,2),s.u0);
shading interp
colorbar
view(2)
axis vis3d
axis tight
