function [util] = smooth4adap(s,u)
%Performs 4th order smoothing od the initial data.
%   s - grid points
%   u - initial condition function handle
%   util - discrete smoothed initial condition corresponding to s


% f4 = ifourier(f4hat(w))
f4 = @(x) (1/36)*(1/2)*...
    (56*x.^3.*sign(x) + (x-3).^3.*(-sign(x-3)) + 12*(x-2).^3.*sign(x-2) -39*(x-1).^3.*sign(x-1) -39*(x+1).^3.*sign(x+1) + 12*(x+2).^3.*sign(x+2) - (x+3).^3.*sign(x+3));
% -1/36 sqrt(?/2) t^3 sgn(t-3)+1/3 sqrt(?/2) t^3 sgn(t-2)-13/12 sqrt(?/2) t^3 sgn(t-1)+7/9 sqrt(2 ?) t^3 sgn(t)-13/12 sqrt(?/2) t^3 sgn(t+1)+1/3 sqrt(?/2) t^3 sgn(t+2)-1/36 sqrt(?/2) t^3 sgn(t+3)+1/4 sqrt(?/2) t^2 sgn(t-3)-sqrt(2 ?) t^2 sgn(t-2)+13/4 sqrt(?/2) t^2 sgn(t-1)-13/4 sqrt(?/2) t^2 sgn(t+1)+sqrt(2 ?) t^2 sgn(t+2)-1/4 sqrt(?/2) t^2 sgn(t+3)-3/4 sqrt(?/2) t sgn(t-3)+2 sqrt(2 ?) t sgn(t-2)-13/4 sqrt(?/2) t sgn(t-1)-13/4 sqrt(?/2) t sgn(t+1)+2 sqrt(2 ?) t sgn(t+2)-3/4 sqrt(?/2) t sgn(t+3)+3/4 sqrt(?/2) sgn(t-3)-4/3 sqrt(2 ?) sgn(t-2)+13/12 sqrt(?/2) sgn(t-1)-13/12 sqrt(?/2) sgn(t+1)+4/3 sqrt(2 ?) sgn(t+2)-3/4 sqrt(?/2) sgn(t+3)
% x = linspace(-3,3,1000);
% h = 1;
% figure(4)
% plot(x,f4(x/h))

% K = 1;
% u = @(x) max(x-K, 0);
% s = linspace(0,4*K,10000);
% h = s(2)-s(1);
% us = u(s);

util = u(s);
for ii = 2:numel(s)-1
    h = mean([abs(s(ii)-s(ii-1)), abs(s(ii+1)-s(ii))]);
    util(ii) = (1/h) * integral(@(x) f4(x/h).*u(s(ii)-x),-3*h,3*h,'AbsTol',1e-32);
end
end
