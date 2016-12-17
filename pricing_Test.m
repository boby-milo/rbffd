clc
clear
close all

dbstop if error

%% Test pricing_init
clear
close all
[~,~,~,~] = pricing_init();

disp('pricing_init: Passed!')

%% Test pricing_grid
clear
close all

% 1D uniform
[contract,parameter,method,grid] = pricing_init();
grid = pricing_grid(grid);

% 2D uniform
[contract,parameter,method,grid] = pricing_init();
grid.dim = 2;
grid = pricing_grid(grid);


%% Test pricing_contract
clear
close all

[contract,~,~,grid] = pricing_init();
grid = pricing_grid(grid);

% 1D
smoothing = 0;
[cnatural] = pricing_contract(contract,grid, smoothing);

smoothing = 1;
[csmooth] = pricing_contract(contract,grid, smoothing);

contract.payoff = 'EUput';
smoothing = 0;
[pnatural] = pricing_contract(contract,grid, smoothing);

smoothing = 1;
[psmooth] = pricing_contract(contract,grid, smoothing);

figure(1)
plot(grid.x,cnatural, grid.x, csmooth, grid.x, pnatural, grid.x, psmooth);

figure(2)
plot(grid.x,cnatural - csmooth, grid.x, pnatural - psmooth);

disp('pricing_contract: Passed!')

%% Test EUcall
clear
close all
%#% Setup
[contract,parameter,method,grid] = pricing_init();

%#% Evaluation
grid = pricing_grid(grid);
grid.u0 = pricing_contract(contract,grid);
grid.W = pricing_differentiation(method,parameter,grid);
[grid.u, ~] = pricing_integration(contract,parameter,method,grid);

figure()
plot(grid.x,grid.u0,grid.x,grid.u)

disp('EUcall: passed!')

%% Test EUbasket
clear
close all
%#% Setup
[contract,parameter,method,grid] = pricing_init();
contract.payoff = 'EUcallBasket';

parameter.sig = [0.3, 0.3];
parameter.rho = 0.5;

grid.dim = 2;
grid.smax = 8*contract.K;

%#% Evaluation
grid = pricing_grid(grid);
grid.u0 = pricing_contract(contract,grid);
grid.W = pricing_differentiation(method,parameter,grid);
% [s.u,m] = pricing_integration(c,p,m,s);
disp('EUcallBasket: passed!')

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
