function [util] = smooth4(s,u,dim)
%Performs 4th order smoothing od the initial data.
%   s - grid points
%   u - initial condition function handle
%   util - discrete smoothed initial condition corresponding to s

if nargin == 2
    dim = 1;
end

% f4 = ifourier(f4hat(w))
f4 = @(x) (1/36)*(1/2)*...
    (56*x.^3.*sign(x) + (x-3).^3.*(-sign(x-3)) + 12*(x-2).^3.*sign(x-2) -39*(x-1).^3.*sign(x-1) -39*(x+1).^3.*sign(x+1) + 12*(x+2).^3.*sign(x+2) - (x+3).^3.*sign(x+3));

switch dim
    case 1
        util = u(s);
        for ii = 2:numel(s)-1
            h = mean([abs(s(ii)-s(ii-1)), abs(s(ii+1)-s(ii))]);
            util(ii) = (1/h) * integral(@(x) f4(x/h).*u(s(ii)-x),-3*h,3*h,'AbsTol',1e-32);
        end
        
    case 2
        util = u(s(:,1),s(:,2));
        h = s(2,2) - s(1,2);
        for ii = 1:numel(s(:,1))
            util(ii) = (1/h^2) * integral2(@(x,y) f4(x/h).*f4(y/h).*u(s(ii,1)-x, s(ii,2)-y),-3*h,3*h,-3*h,3*h);
        end
end
end
