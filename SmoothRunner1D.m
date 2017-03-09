clc
clear
% close all

dbstop if error

Kmul = 4;


runnumber=1;

N = [40, 80, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2000, 2400, 2800, 3600, 4000, 8000, 40000];
M = 5*N;

Kx = 1;

legnames = {'FD2half','FD4half','FD2',...
    'FD4','FD4smooth','RBFFDadap',...
    'RBFFDadap-smooth','RBFFDreg-phs',...
    'RBFFDreg-smooth-phs','RBFFDadap-phs',...
    'RBFFDadap-smooth-phs'};

%% Run
for ii = 1:numel(N)
    
    g=15;
    n=5;
    
    %FD
    F1 = parfeval(@BSeuCall1D_FD2half,8,N(ii),M(ii),Kmul);
    F2 = parfeval(@BSeuCall1D_FD4half,8,N(ii),M(ii),Kmul);
    F3 = parfeval(@BSeuCall1D_FD2,8,N(ii),M(ii),Kmul);
    F4 = parfeval(@BSeuCall1D_FD4,8,N(ii),M(ii),Kmul);
    F11 = parfeval(@BSeuCall1D_FD4smooth,8,N(ii),M(ii),Kmul);
    
    %RBF
    F5 = parfeval(@BSeuCall1D_RBFFDadap,8,N(ii),g,n,M(ii),Kmul);
    F6 = parfeval(@BSeuCall1D_RBFFDadap_smooth,8,N(ii),g,n,M(ii),Kmul);
    
    %PHS
    d = 5;
    p = 5;
    nm = 5;
    
    F7 = parfeval(@BSeuCall1D_RBFFDreg_phs,8,N(ii),p,d,nm,M(ii),Kmul);
    F8 = parfeval(@BSeuCall1D_RBFFDreg_smooth_phs,8,N(ii),p,d,nm,M(ii),Kmul);
    F9 = parfeval(@BSeuCall1D_RBFFDadap_phs,8,N(ii),g,p,d,nm,M(ii),Kmul);
    F10 = parfeval(@BSeuCall1D_RBFFDadap_smooth_phs,8,N(ii),g,p,d,nm,M(ii),Kmul);
    
    
    [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),~,~,~] = fetchOutputs(F1);
    disp([int2str(ii),' F1 ', num2str(timF1)])
    [uF2{ii},errF2{ii},timF2,xF2{ii},dxF2(ii),~,~,~] = fetchOutputs(F2);
    disp([int2str(ii),' F2 ', num2str(timF2)])
    [uF3{ii},errF3{ii},timF3,xF3{ii},dxF3(ii),~,~,~] = fetchOutputs(F3);
    disp([int2str(ii),' F3 ', num2str(timF3)])
    
    [uF4{ii},errF4{ii},timF4,xF4{ii},dxF4(ii),~,~,~] = fetchOutputs(F4);
    disp([int2str(ii),' F4 ', num2str(timF4)])
    
    [uF5{ii},errF5{ii},timF5,xF5{ii},dxF5(ii),~,~,~] = fetchOutputs(F5);
    disp([int2str(ii),' F5 ', num2str(timF5)])
    
    [uF6{ii},errF6{ii},timF6,xF6{ii},dxF6(ii),~,~,~] = fetchOutputs(F6);
    disp([int2str(ii),' F6 ', num2str(timF6)])
    
    [uF7{ii},errF7{ii},timF7,xF7{ii},dxF7(ii),~,~,~] = fetchOutputs(F7);
    disp([int2str(ii),' F7 ', num2str(timF7)])
    
    [uF8{ii},errF8{ii},timF8,xF8{ii},dxF8(ii),~,~,~] = fetchOutputs(F8);
    disp([int2str(ii),' F8 ', num2str(timF8)])
    
    [uF9{ii},errF9{ii},timF9,xF9{ii},dxF9(ii),~,~,~] = fetchOutputs(F9);
    disp([int2str(ii),' F9 ', num2str(timF9)])
    
    [uF10{ii},errF10{ii},timF10,xF10{ii},dxF10(ii),~,~,~] = fetchOutputs(F10);
    disp([int2str(ii),' F10 ', num2str(timF10)])
    
     [uF11{ii},errF11{ii},timF11,xF11{ii},dxF11(ii),~,~,~] = fetchOutputs(F11);
    disp([int2str(ii),' F11 ', num2str(timF11)])
    
    
    indreg=[];
    for jj = 1:length(xF1{ii})
        %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
        if xF1{ii}(jj) >= 1/3*Kx && xF1{ii}(jj) <= 5/3*Kx
            indreg=[indreg jj];
        end
    end
    
    xind = xF1{ii}(indreg);
    
    
    
    errF1crop{ii} = errF1{ii}(indreg);
    errF2crop{ii} = errF2{ii}(indreg);
    errF3crop{ii} = errF3{ii}(indreg);
    errF4crop{ii} = errF4{ii}(indreg);
    errF5crop{ii} = errF5{ii}(indreg);
    errF6crop{ii} = errF6{ii}(indreg);
    errF7crop{ii} = errF7{ii}(indreg);
    errF8crop{ii} = errF8{ii}(indreg);
    errF9crop{ii} = errF9{ii}(indreg);
    errF10crop{ii} = errF10{ii}(indreg);
    errF11crop{ii} = errF11{ii}(indreg);
    
    
    
    erF1(ii) = max(abs(errF1crop{ii}));
    erF2(ii) = max(abs(errF2crop{ii}));
    erF3(ii) = max(abs(errF3crop{ii}));
    erF4(ii) = max(abs(errF4crop{ii}));
    erF5(ii) = max(abs(errF5crop{ii}));
    erF6(ii) = max(abs(errF6crop{ii}));
    erF7(ii) = max(abs(errF7crop{ii}));
    erF8(ii) = max(abs(errF8crop{ii}));
    erF9(ii) = max(abs(errF9crop{ii}));
    erF10(ii) = max(abs(errF10crop{ii}));
    erF11(ii) = max(abs(errF11crop{ii}));
    
    figure(runnumber)
    clf
    loglog(dxF1,erF1,'--','LineWidth',2)
    hold on
    loglog(dxF2,erF2,'--','LineWidth',2)
    loglog(dxF3,erF3,'-','LineWidth',2)
    loglog(dxF4,erF4,'-','LineWidth',2)
    loglog(dxF11,erF11,'-.','LineWidth',2)
    
    loglog(dxF5,erF5,'-sq','MarkerSize',12,'LineWidth',2)
    loglog(dxF6,erF6,'-.sq','MarkerSize',12,'LineWidth',2)
    
    loglog(dxF7,erF7,'-^','MarkerSize',12,'LineWidth',2)
    loglog(dxF8,erF8,'-.*','MarkerSize',12,'LineWidth',2)
    loglog(dxF9,erF9,'-^','MarkerSize',12,'LineWidth',2)
    loglog(dxF10,erF10,'-.*','MarkerSize',12,'LineWidth',2)
    legend(legnames);
    %     title(['M=',num2str(M)]);
        xlim([1e-4, 1e-1]);
            ylim([1e-11, 1e-1]);
    title('BSeuCall1D')
    xlabel('h')
    ylabel('\Deltau_{max}')
    set(gca,'fontsize',18)
    drawnow
end