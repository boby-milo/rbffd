function [s,N,indres,indin,indcf,indff,inddy0] = HSTeuCall_grid_uniform(Nx, Kx, sres, slim, pl)

x_min = slim(1);
x_max = slim(2);
y_min = slim(3);
y_max = slim(4);

xres = sres(:,1); yres = sres(:,2);

dx = 1/Nx;

%% Uniform
xc = Kx;
yc = yres(1);

xdo = linspace(x_min,x_max,Nx);
ydo = zeros(size(xdo)) + y_min;
xup = linspace(x_min,x_max,Nx);
yup = zeros(size(xup)) + y_max;
yle = linspace(y_min,y_max,Nx);
xle = zeros(size(yle)) + x_min;
yri = linspace(y_min,y_max,Nx);
xri = zeros(size(yri)) + x_max;

a = 0.75; b = 0.25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius = @(x,y) dx * ( (x-xc).^2/a^2 + (y-yc).^2/b^2 + 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box = [x_min-dx,x_max+dx,y_min-dx,y_max+dx];

corners = [x_min y_min; x_max y_min; x_max y_max; x_min y_max];
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
indcf = find(xvec==0)';
indff = find(xvec==x_max)'; %, Ny:Ny:N-Ny];
indin = ind;
indup = find(yvec==y_max)';
inddo = find(yvec==y_min)';

indyup = indup(2:end-1);
indydo = inddo(2:end-1);
inddy0 = [];
% inddy0 = [indyup, indydo];
% inddy0 = indyup;

indin([indff,indcf,inddy0]) = [];
% indin([indff,indcf]) = [];

indres = [find(xvec == xres(1) & yvec == yres(1)),...
    find(xvec == xres(2) & yvec == yres(2)),...
    find(xvec == xres(3) & yvec == yres(3))];

% % equidistant
% x = transpose(linspace(x_min,x_max,Nx));
% Ny = round(Nx*((y_max-y_min)/x_max));
% y = transpose(linspace(y_min,y_max,Ny));
%
% x = [x;xress];
% y = [y;yress];
%
% [X,Y] = meshgrid(x,y);
%
% xvec = X(:);
% yvec = Y(:);

% % halton
% p1 = haltonset(1);
% p2 = haltonset(2);
%
% xdo = net(p1,Nx);
% ydo = zeros(size(xdo))+y_min;
% xup = net(p1,Nx);
% yup = ones(size(xup));
% yle = net(p1,Nx);
% xle = zeros(size(yle));
% yri = net(p1,Nx);
% xri = ones(size(yri));
%
% s = net(p2,Nx^2);
% xvec = [s(:,1);xdo;xup;xle;xri];
% yvec = [s(:,2);ydo;yup;yle;yri];

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
