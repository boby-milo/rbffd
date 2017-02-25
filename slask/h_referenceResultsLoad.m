clc
clear
close all


warning off
dbstop if error

leg{1} = 'BSeuCall2Dbasket_FD.mat';
leg{2} = 'BSeuCall2Dbasket_RBFFDreg_phs.mat';

load(leg{1})

for jj = 1:numel(leg)
    load(leg{jj})
    
    for ii = 1:numel(N)
        errabs(ii) = max(abs(err{ii}));
        dxabs(ii) = dx{ii};
    end
    
    figure(1)
    loglog(dxabs,errabs,'-*')
    hold on
end

xo = 10.^(-2:0.1:0);
yo2 = (1/10)*xo.^2;
yo4 = (1/10)*xo.^4;
loglog(xo,yo2,'k--')
loglog(xo,yo4,'k--')
xlim([1e-2,1e0]);
ylim([1e-9,1e-1]);
legend(leg)