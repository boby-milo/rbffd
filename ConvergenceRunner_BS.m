clc
clear
% close all
dbstop if error
warning off

cd('./BS2D')

i = linspace(0.75*1.699,1.2*2.3010,25);
Nx = round(10.^i);

u1 = nan(1,numel(Nx));
u2 = u1;
u3 = u1;

% load(fileName);

pl = 1;
for ii = 1:numel(Nx)
    ii
    p = gcp();
    
    [ures1,errres1,tim1,N1] = BSamPut2D_FD2(Nx(ii))
    t1(ii) = tim1;
    err1(ii) = max(abs(errres1));
    dx1(ii) = 1/sqrt(N1);
    
%     [ures1,errres1,tim1,cond1,N1] = BSeuCall2D_FD2(Nx(ii))
%     cnd1(ii) = cond1;
%     t1(ii) = tim1;
%     err1(ii) = max(abs(errres1));
%     dx1(ii) = 1/sqrt(N1);
    
%     [ures1,errres1,tim1,cond1,N1] = BSeuCall2D_RBFFD_cartesian(Nx(ii),pl)
%     cnd1(ii) = cond1;
%     t1(ii) = tim1;
%     err1(ii) = max(abs(errres1));
%     dx1(ii) = 1/sqrt(N1);
%     
%     [ures2,errres2,tim2,cond2,N2] = BSeuCall2D_RBFFD_adap(Nx(ii),pl)
%     cnd2(ii) = cond2;
%     t2(ii) = tim2;
%     err2(ii) = max(abs(errres2));
%     dx2(ii) = 1/sqrt(N2);
    
    [ures3,errres3,tim3,cond3,N3] = BSamPut2D_RBFFD(2*Nx(ii),pl)
    cnd3(ii) = cond3;
    t3(ii) = tim3;
    err3(ii) = max(abs(errres3));
    dx3(ii) = 1/sqrt(N3);
    
    
    figure(101)
    clf
    loglog(dx1,err1,dx3,err3) 
    hold on
    legend('cartesian', 'adapted', 'uniform')
    xlabel('\Deltax')
    ylabel('error')
    drawnow
    
    figure(102)
    clf
    loglog(t1,err1,t3,err3) 
    hold on
    legend('cartesian', 'adapted', 'uniform')
    xlabel('time')
    ylabel('error')
    drawnow
    
%     figure(103)
%     clf
%     loglog(dx1,cnd1) 
%     hold on
%     legend('cartesian', 'adapted', 'uniform')
%     xlabel('\Deltax')
%     ylabel('\kappa(A)')
%     drawnow
    
end

cd('../')