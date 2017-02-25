clc
clear
close all


warning off
dbstop if error

% ppar = gcp();
% pctRunOnAll warning off

h_referenceResultsLoad

Kmul = 8;

N = [40, 80, 120, 160, 200, 240, 280, 320];
M = 4*N;

xo = 10.^(-2:0.1:0);
yo2 = (1/10)*xo.^2;
yo4 = (1/10)*xo.^4;

for ii = 1:numel(N)
    ii
    
%     F0 = parfeval(@BSeuCallbasket2D_FD,7,N(ii),M(ii),Kmul);
   
    d = 5;
    p = 7;
    g = 10;
    nm = 2.5;
    
    F0 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M(ii),Kmul,nm);
    
   
%     F1 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M(ii),Kmul,nm);
%     F3 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs_smooth,8,N(ii),p,d,M(ii),Kmul,nm);  
    
    [u{ii},err{ii},tim{ii},x{ii},dx{ii},Ntot{ii},~] = fetchOutputs(F0);
    
    
    %     disp([int2str(ii),' F0 ', num2str(timF0)])
    %     disp([NtotF0(ii),d,p]);
    er(ii) = norm(err{ii},Inf);
    dxplot(ii) = dx{ii};
    %     disp(erF0(ii));
%     
%     [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),nF1(ii),NtotF1(ii),~] = fetchOutputs(F1);
%     disp([int2str(ii),' F1 ', num2str(timF1)])
%     disp([NtotF1(ii),d,p,nF1(ii)]);
%     erF1(ii) = norm(errF1{ii},Inf);
%     disp(erF1(ii));
    
    f=figure(1);
%     clf
    
    loglog(dxplot,er,'-sq')
    
    
%     loglog(xo,yo2,'k--')
%     loglog(xo,yo4,'k--')
    xlim([1e-2,1e0]);
    ylim([1e-9,1e-1]);
    
    
%     title(runtitle)
%     xlabel('h')
%     ylabel('\Deltau_{max}')
%     set(gca,'fontsize',18)
    drawnow
end

save('BSeuCall2Dbasket_RBFFDreg_phs_try.mat', 'u', 'err', 'tim', 'x', 'dx', 'Ntot', 'N', 'M', 'Kmul')