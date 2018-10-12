function [bdyr] = discretize_bdy (corners,radius,close_p)

% Output parameter
%   bdyr        Array with 3 columns containing x,y,r for all bdy points
%
% Input parameters
%   corners     List of x,y coordinates for all corners; 2 columns.
%   radius      Handle for function returning radii
%   close_p     true:  code closes polygon
%               false: polygon is left open

% The routine will add extra boundary markers whenever a side of the 
% polygonal boundary is longer than about 1.2 times the local grain radius.

z_init = complex(corners(:,1),corners(:,2));
                % Algorithm becomes somewhat simpler if using complex
                % arithmetic.
l_z_init = length(z_init);      % Nr of initial polygon corner points

dotmax = 1e+4;                  % Max allowed number of final dots
z      = zeros(dotmax,1);       % Vector to store boundary nodes
m      = 0;                     % Counter in z-vector

z_init(l_z_init+1) = z_init(1); % Repeat first element
for k = 1:l_z_init              % Loop over the sides of the initial polygon
    dir = sign(z_init(k+1)-z_init(k));   % Unit number pointing along side
    m  = m+1;
    mk = m;                     % Remember start point along present side
    z(m) = z_init(k);           % Store boundary nodes in z
    if k==l_z_init & ~close_p; break;end
    r  = radius(real(z(m)),imag(z(m)));
    while abs(z(m)-z_init(k+1)) > 1.2*r
        m = m+1;
        z(m) = z(m-1)+r*dir;
        r  = radius(real(z(m)),imag(z(m)));
    end
    
    % Adjust nodes along the present side to match nicely at far end
    np = m - mk;               % Number of points to adjust
    if np > 0
        z_mm = z(m)+r*dir;   % This would have been the next point
        z(mk+1:m) = z(mk)+(z(mk+1:m)-z(mk))*... % Ajust along side
            abs((z_init(k+1)-z(mk))/(z_mm-z(mk)));
    end
end
dotnr = m;

bdyr    = [real(z(1:m)),imag(z(1:m)),radius(real(z(1:m)),imag(z(1:m)))];