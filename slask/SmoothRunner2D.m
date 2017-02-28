clc
clear
% close all

warning off
dbstop if error

ppar = gcp();
pctRunOnAll warning off

Kmul = 8;

runnumber = 1;

N = [40, 80, 120, 160, 200, 240, 280, 320];
M = 4*N;

xo = 10.^(-2:0.1:0);
yo2 = (1/10)*xo.^2;
yo4 = (1/10)*xo.^4;

for ii = 1:numel(N)
    
    
    F0 = parfeval(@BSeuCallbasket2D_FD,7,N(ii),M(ii),Kmul);
    
    nm=2.2;
    
    g = 5;
    
    d = 3;
    p = 5;
    
    F1 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M(ii),Kmul,nm);
    F3 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs_smooth,8,N(ii),p,d,M(ii),Kmul,nm);
    F5 = parfeval(@BSeuCall2Dbasket_RBFFDadap_phs,8,N(ii),p,d,M(ii),Kmul,g,nm);
    F7 = parfeval(@BSeuCall2Dbasket_RBFFDadap_phs_smooth,8,N(ii),p,d,M(ii),Kmul,g,nm);
    
    d = 5;
    p = 7;
    F2 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M(ii),Kmul,nm);
    F4 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs_smooth,8,N(ii),p,d,M(ii),Kmul,nm);
    F6 = parfeval(@BSeuCall2Dbasket_RBFFDadap_phs,8,N(ii),p,d,M(ii),Kmul,g,nm);
    F8 = parfeval(@BSeuCall2Dbasket_RBFFDadap_phs_smooth,8,N(ii),p,d,M(ii),Kmul,g,nm);
    
    
    [uF0{ii},errF0{ii},timF0,xF0{ii},dxF0(ii),NtotF0(ii),~] = fetchOutputs(F0);
    %     disp([int2str(ii),' F0 ', num2str(timF0)])
    %     disp([NtotF0(ii),d,p]);
    erF0(ii) = norm(errF0{ii},Inf);
    %     disp(erF0(ii));
    
    [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),nF1(ii),NtotF1(ii),~] = fetchOutputs(F1);
    disp([int2str(ii),' F1 ', num2str(timF1)])
    disp([NtotF1(ii),d,p,nF1(ii)]);
    erF1(ii) = norm(errF1{ii},Inf);
    disp(erF1(ii));
    
    [uF2{ii},errF2{ii},timF2,xF2{ii},dxF2(ii),nF2(ii),NtotF2(ii),~] = fetchOutputs(F2);
    disp([int2str(ii),' F2 ', num2str(timF2)])
    disp([NtotF2(ii),d,p,nF2(ii)]);
    erF2(ii) = norm(errF2{ii},Inf);
    disp(erF2(ii));
    
    [uF3{ii},errF3{ii},timF3,xF3{ii},dxF3(ii),nF3(ii),NtotF3(ii),~] = fetchOutputs(F3);
    disp([int2str(ii),' F3 ', num2str(timF3)])
    disp([NtotF3(ii),d,p,nF3(ii)]);
    erF3(ii) = norm(errF3{ii},Inf);
    disp(erF3(ii));
    
    [uF4{ii},errF4{ii},timF4,xF4{ii},dxF4(ii),nF4(ii),NtotF4(ii),~] = fetchOutputs(F4);
    disp([int2str(ii),' F4 ', num2str(timF4)])
    disp([NtotF4(ii),d,p,nF4(ii)]);
    erF4(ii) = norm(errF4{ii},Inf);
    disp(erF4(ii));
    
    [uF5{ii},errF5{ii},timF5,xF5{ii},dxF5(ii),nF5(ii),NtotF5(ii),~] = fetchOutputs(F5);
    disp([int2str(ii),' F5 ', num2str(timF5)])
    disp([NtotF5(ii),d,p,nF5(ii)]);
    erF5(ii) = norm(errF5{ii},Inf);
    disp(erF5(ii));
    
    [uF6{ii},errF6{ii},timF6,xF6{ii},dxF6(ii),nF6(ii),NtotF6(ii),~] = fetchOutputs(F6);
    disp([int2str(ii),' F6 ', num2str(timF6)])
    disp([NtotF6(ii),d,p,nF6(ii)]);
    erF6(ii) = norm(errF6{ii},Inf);
    disp(erF6(ii));
    
    [uF7{ii},errF7{ii},timF7,xF7{ii},dxF7(ii),nF7(ii),NtotF7(ii),~] = fetchOutputs(F7);
    disp([int2str(ii),' F7 ', num2str(timF7)])
    disp([NtotF7(ii),d,p, nF7(ii)]);
    erF7(ii) = norm(errF7{ii},Inf);
    disp(erF7(ii));
    
    [uF8{ii},errF8{ii},timF8,xF8{ii},dxF8(ii),nF8(ii),NtotF8(ii),~] = fetchOutputs(F8);
    disp([int2str(ii),' F8 ', num2str(timF8)])
    disp([NtotF8(ii),d,p, nF8(ii)]);
    erF8(ii) = norm(errF8{ii},Inf);
    disp(erF8(ii));
    
    f=figure(runnumber);
    clf
    loglog(dxF0,erF0,'-sq')
    hold on
    loglog(dxF1,erF1,'-^')
    loglog(dxF2,erF2,'-^')
    loglog(dxF3,erF3,'-o')
    loglog(dxF4,erF4,'-o')
    loglog(dxF5,erF5,'-*')
    loglog(dxF6,erF6,'-*')
        loglog(dxF7,erF7,'-')
        loglog(dxF8,erF8,'-')
    
    loglog(xo,yo2,'k--')
    loglog(xo,yo4,'k--')
    xlim([1e-2,1e0]);
    ylim([1e-9, 1e-1]);
    
    title('BSeuCall2Dbasket')
    legend('FD2\_reg\_square',... %F0
        'PHSd3p5\_reg\_triangle',... %F1
        'PHSd5p7\_reg\_triangle',... %F2
        'PHSd3p5\_reg\_triangle\_smooth',... %F3
        'PHSd5p7\_reg\_triangle\_smooth',... %F4
        'PHSd3p5\_adap\_triangle',... %F5
        'PHSd5p7\_adap\_triangle',... %F6
        'PHSd3p5\_adap\_triangle\_smooth',... %F7
        'PHSd5p7\_adap\_triangle\_smooth') %F8
    xlabel('h')
    ylabel('\Deltau_{max}')
    set(gca,'fontsize',18)
    drawnow
end