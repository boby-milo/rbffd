function [s,N,indres,indin,indcf,indff] = BSeuCall2D_grid_uniform(Nx, Kx, sres, slim, pl)

x_min = slim(1);
x_max = slim(2);
y_min = slim(3);
y_max = slim(4);

xres = sres(:,1); yres = sres(:,2);

dx = 1/Nx;

%% Uniform
xc = Kx;
yc = yres(1);

% xdo = linspace(x_min,x_max,Nx);
% ydo = zeros(size(xdo)) + y_min;
% xup = linspace(x_min,x_max,Nx);
% yup = zeros(size(xup)) + y_max;
% yle = linspace(y_min,y_max,Nx);
% xle = zeros(size(yle)) + x_min;
% yri = linspace(y_min,y_max,Nx);
% xri = zeros(size(yri)) + x_max;

a = 0.25; b = 0.75; d = pi/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius = @(x,y) dx * (   ( (x-xc)*cos(d) + (y-yc)*sin(d) ).^2 / a^2   +   ( (x-xc)*sin(d) - (y-yc)*cos(d) ).^2 / b^2 + 1   );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box = [x_min-dx,x_max+dx,y_min-dx,y_max+dx];

corners = [x_min y_min; x_max y_min; x_min y_max];
boundaries = discretize_bdy(corners,radius,true);
boundaries = [boundaries(:,1:2); [xres yres]];

s = node_drop(box, radius);
s = repel(s, boundaries, corners, @(x,y) radius(x,y));
s = [s; boundaries];
xvec = s(:,1);
yvec = s(:,2);

ind = find(xvec<=x_max & yvec<=y_max & xvec>=x_min & yvec>=y_min);
xvec = xvec(ind);
yvec = yvec(ind);

s = unique([xvec, yvec], 'rows');
xvec = s(:,1);
yvec = s(:,2);

N = numel(xvec);
Nx_res = sqrt(N);

ind = 1:N;
indcf = find(xvec==0 & yvec ==0)';
indff = find(abs(yvec + xvec - x_max) <= 1e-12)'; %, Ny:Ny:N-Ny];
indin = ind;
% indup = find(yvec==y_max)';
% inddo = find(yvec==y_min)';

% indyup = indup(2:end-1);
% indydo = inddo(2:end-1);
% inddy0 = [];
% inddy0 = [indyup, indydo];
% inddy0 = indyup;

indin([indff,indcf]) = [];
% indin([indff,indcf]) = [];

indres = [find(xvec == xres(1) & yvec == yres(1)),...
    find(xvec == xres(2) & yvec == yres(2)),...
    find(xvec == xres(3) & yvec == yres(3))];

if pl
    figure(1)
    clf
    plot(xvec,yvec,'k.')
    axis equal
    axis tight
    hold on
%     plot(xvec(indin),yvec(indin),'ko')
    plot(xvec(indcf),yvec(indcf),'^')
    plot(xvec(indff),yvec(indff),'sq')
    plot(xvec(indres),yvec(indres),'*')
    drawnow
end
end
