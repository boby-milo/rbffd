function [util] = smooth4adap(s,dgs,u,hloc)

% f4 = ifourier(f4hat(w))
f4 = @(x) (1/36)*(1/2)*...
    (56*x.^3.*sign(x) + (x-3).^3.*(-sign(x-3)) + 12*(x-2).^3.*sign(x-2) -39*(x-1).^3.*sign(x-1) -39*(x+1).^3.*sign(x+1) + 12*(x+2).^3.*sign(x+2) - (x+3).^3.*sign(x+3));

util = u(s(:,1),s(:,2));
h = min(hloc);
% indreg = [];
for ii = 1:numel(dgs)
    udgs = (1/h^2) * integral2(@(x,y) f4(x/h).*f4(y/h).*u(s(dgs{ii}(1),1)-x, s(dgs{ii}(1),2)-y),-3*h,3*h,-3*h,3*h,'AbsTol',1e-12);
    util(dgs{ii}) = udgs;
%     indreg = [indreg; dgs{ii}];
end
end
