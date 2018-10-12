clc
clear
close

cd('BSeuCall2D')

q = [1,3,5,7,9];

r = linspace(0,1,200);

for ii = 1:numel(q)
    phi = phs(q(ii),r,'0',1)
    
    figure(1)
    plot(r,phi)
    hold on
    axis equal
    axis tight
end

