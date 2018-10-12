clc
clear
% close all

warning off

cd('./HSTeuCall')

fileName = 'HSTeuCall_macbook';
i = linspace(1.75,3,25);
Nx = round(10.^i);

u1 = nan(1,numel(Nx));
u2 = u1;
u3 = u1;

% load(fileName);

pl = 1;
for ii = 1:numel(Nx)
    ii
    p = gcp();
    
    [ures1,errres1,tim1,cond1,N1] = HSTeuCall_RBFFD_cartesian(round(Nx(ii)*0.4),0)
    cnd1(ii) = cond1;
    t1(ii) = tim1;
    err1(ii) = max(abs(errres1));
    dx1(ii) = 1/sqrt(N1);
    
    [ures2,errres2,tim2,cond2,N2] = HSTeuCall_RBFFD_adap(round(Nx(ii)*0.4),0)
    cnd2(ii) = cond2;
    t2(ii) = tim2;
    err2(ii) = max(abs(errres2));
    dx2(ii) = 1/sqrt(N2);
    
    [ures3,errres3,tim3,cond3,N3] = HSTeuCall_RBFFD(Nx(ii),pl)
    cnd3(ii) = cond3;
    t3(ii) = tim3;
    err3(ii) = max(abs(errres3));
    dx3(ii) = 1/sqrt(N3);
    
    
    figure(201)
    clf
    loglog(dx1,err1,dx2,err2,dx3,err3)
    grid on
    hold on
    legend('cartesian', 'adapted', 'uniform')
    xlabel('\Deltax')
    ylabel('error')
    ylim([1e-4,1])
    drawnow
    
    figure(202)
    clf
    loglog(t1,err1,t2,err2,t3,err3)
    grid on
    hold on
    legend('cartesian', 'adapted', 'uniform')
    xlabel('time')
    ylabel('error')
    drawnow
    
    figure(203)
    clf
    loglog(dx1,cnd1,dx2,cnd2,dx3,cnd3) 
    grid on
    hold on
    legend('cartesian', 'adapted', 'uniform')
    xlabel('\Deltax')
    ylabel('\kappa(A)')
    drawnow
    
end

cd('../')