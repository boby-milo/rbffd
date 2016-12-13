clc
clear
close all

dbstop if error

%% Test default_init
clear
%% Test make_grid
clear
%% Test payoff_function
clear
close all

[contract,~,~,grid] = default_init();
grid = make_grid(grid);

% 1D
smoothing = 0;
[unatural] = payoff_function(contract,grid, smoothing);

smoothing = 1;
[usmooth] = payoff_function(contract,grid, smoothing);

figure(1)
plot(grid.x,unatural, grid.x, usmooth);

figure(2)
plot(grid.x,unatural-usmooth);

disp('payoff_function: Passed!')

%% Test EUcall
clear
close all
%#% Setup
[contract,parameter,method,grid] = default_init();

%#% Evaluation
grid = make_grid(grid);
grid.u0 = payoff_function(contract,grid);
grid.W = differentiation_matrix(method,parameter,grid);
[grid.u, ~] = time_integrate(contract,parameter,method,grid);

figure()
plot(grid.x,grid.u0,grid.x,grid.u)

disp('EUcall test passed!')

%% Test EUbasket
clear
close all
%#% Setup
[contract,parameter,method,grid] = default_init();
contract.payoff = 'EUcallBasket';

parameter.sig = [0.3, 0.3];
parameter.rho = 0.5;

grid.dim = 2;
grid.smax = 8*contract.K;

%#% Evaluation
grid = make_grid(grid);
grid.u0 = payoff_function(contract,grid);
grid.W = differentiation_matrix(method,parameter,grid);
% [s.u,m] = time_integrate(c,p,m,s);
disp('EUcallBasket test passed!')

figure()
clf
plot(grid.x(grid.indff,1),grid.x(grid.indff,2),'b^','MarkerFaceColor','auto')
hold on
plot(grid.x(grid.indcf,1),grid.x(grid.indcf,2),'rsq','MarkerFaceColor','auto')
plot(grid.x(grid.indin,1),grid.x(grid.indin,2),'ko','MarkerFaceColor','auto')
plot(grid.x(:,1),grid.x(:,2),'k.')
axis equal
axis tight
hold off

figure()
tri = delaunay(grid.x(:,1),grid.x(:,2));
trisurf(tri, grid.x(:,1),grid.x(:,2),grid.u0);
shading interp
colorbar
view(2)
axis vis3d
axis tight
