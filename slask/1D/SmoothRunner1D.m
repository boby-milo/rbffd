clc
clear
% close all

dbstop if error

Kmul = 4;
M = 1000; 

runnumber=1;

N = [25, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600];

Kx = 1;

%% Run
for ii = 1:numel(N)
    
    n=5;
    F0 = parfeval(@BSeuCall1D_RBFFDadapsmooth,7,N(ii),n,M,Kmul);
    
    F1 = parfeval(@BSeuCall1D_FD2half,8,N(ii),M,Kmul);
    F2 = parfeval(@BSeuCall1D_FD4half,8,N(ii),M,Kmul);
    F3 = parfeval(@BSeuCall1D_FD2,7,N(ii),M,Kmul);
    
    
    
    
    d = 5;
    p = 5; 
%     F4 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
    F4 = parfeval(@BSeuCall1D_RBFFDreg_phs_smooth,7,N(ii),p,d,M,Kmul);
    
%     F5 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
    F5 = parfeval(@BSeuCall1D_RBFFDreg_phs_adap_smooth,7,N(ii),p,d,M,Kmul);
    
    p = 7; 
%     F6 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
    F6 = parfeval(@BSeuCall1D_RBFFDreg_phs_adap_smooth,7,N(ii),p,d,M,Kmul);
    
    [uF0{ii},errF0{ii},timF0,xF0{ii},dxF0(ii),~,~] = fetchOutputs(F0);
    disp([int2str(ii),' F0 ', num2str(timF0)])
    [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),~,~,~] = fetchOutputs(F1);
    disp([int2str(ii),' F1 ', num2str(timF1)])
    [uF2{ii},errF2{ii},timF2,xF2{ii},dxF2(ii),~,~,~] = fetchOutputs(F2);
    disp([int2str(ii),' F2 ', num2str(timF2)])
    [uF3{ii},errF3{ii},timF3,xF3{ii},dxF3(ii),~,~] = fetchOutputs(F3);
    disp([int2str(ii),' F3 ', num2str(timF3)])
    
    [uF4{ii},errF4{ii},timF4,xF4{ii},dxF4(ii),~,~] = fetchOutputs(F4);
    disp([int2str(ii),' F4 ', num2str(timF4)])
    
    [uF5{ii},errF5{ii},timF5,xF5{ii},dxF5(ii),~,~] = fetchOutputs(F5);
    disp([int2str(ii),' F5 ', num2str(timF5)])
    
    [uF6{ii},errF6{ii},timF6,xF6{ii},dxF6(ii),~,~] = fetchOutputs(F6);
    disp([int2str(ii),' F6 ', num2str(timF6)])
    
    
    indreg=[];
    for jj = 1:length(xF1{ii})
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xF1{ii}(jj) >= 1/3*Kx && xF1{ii}(jj) <= 5/3*Kx
            indreg=[indreg jj];
        end
    end
    
    xind = xF1{ii}(indreg);
    
    errF5crop{ii} = errF5{ii}(indreg);
    errF6crop{ii} = errF6{ii}(indreg);
    errF0crop{ii} = errF0{ii}(indreg);
    errF4crop{ii} = errF4{ii}(indreg);
    errF1crop{ii} = errF1{ii}(indreg);
    errF2crop{ii} = errF2{ii}(indreg);
    errF3crop{ii} = errF3{ii}(indreg);
    
    erF5(ii) = max(abs(errF5crop{ii}));
    erF6(ii) = max(abs(errF6crop{ii}));
    erF0(ii) = max(abs(errF0crop{ii}));
    erF4(ii) = max(abs(errF4crop{ii}));
    erF1(ii) = max(abs(errF1crop{ii}));
    erF2(ii) = max(abs(errF2crop{ii}));
    erF3(ii) = max(abs(errF3crop{ii}));
    
    figure(runnumber)
    clf
    loglog(dxF0,erF0,'-*')
    hold on
    loglog(dxF1,erF1,'-*')
    loglog(dxF2,erF2,'-*')
    loglog(dxF3,erF3,'-*')
    loglog(dxF4,erF4,'-*')
    loglog(dxF5,erF5,'-*')
    loglog(dxF6,erF6,'-*')
    legend('F0', 'F1', 'F2', 'F3', 'F4','F5', 'F6');
    title(['M=',num2str(M)]);
    %     xlim([1e-4, 1e-1]);
%         ylim([1e-9, 1e-1]);
    drawnow
end