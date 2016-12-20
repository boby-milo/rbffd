clc
clear
% close all

warning off
dbstop if error

ppar = gcp();
pctRunOnAll warning off

Kmul = 8;
M = 100;

runnumber = 1;

N = [20, 40, 60, 80, 100, 120, 140, 160];


for ii = 1:numel(N)
    
    
    F0 = parfeval(@BSeuCallbasket2D_FD,7,N(ii),M,Kmul);
    
    d = 5;
    p = 5;
    F1 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M,Kmul);
    F3 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs_smooth,8,N(ii),p,d,M,Kmul);
    
    d = 7;
    p = 7;
    F2 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs,8,N(ii),p,d,M,Kmul);
    F4 = parfeval(@BSeuCall2Dbasket_RBFFDreg_phs_smooth,8,N(ii),p,d,M,Kmul);
    
    
    
    [uF0{ii},errF0{ii},timF0,xF0{ii},dxF0(ii),NtotF0(ii),~] = fetchOutputs(F0);
%     disp([int2str(ii),' F0 ', num2str(timF0)])
%     disp([NtotF0(ii),d,p]);
    erF0(ii) = norm(errF0{ii},Inf);
%     disp(erF0(ii));
    
    [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),nF1(ii),NtotF1(ii),~] = fetchOutputs(F1);
    disp([int2str(ii),' F1 ', num2str(timF1)])
    disp([NtotF1(ii),d,p]);
    erF1(ii) = norm(errF1{ii},Inf);
    disp(erF1(ii));
    
    [uF2{ii},errF2{ii},timF2,xF2{ii},dxF2(ii),nF2(ii),NtotF2(ii),~] = fetchOutputs(F2);
    disp([int2str(ii),' F2 ', num2str(timF2)])
    disp([NtotF2(ii),d,p]);
    erF2(ii) = norm(errF2{ii},Inf);
    disp(erF2(ii));
    
    [uF3{ii},errF3{ii},timF3,xF3{ii},dxF3(ii),nF3(ii),NtotF3(ii),~] = fetchOutputs(F3);
    disp([int2str(ii),' F3 ', num2str(timF3)])
    disp([NtotF3(ii),d,p]);
    erF3(ii) = norm(errF3{ii},Inf);
    disp(erF3(ii));
    
    [uF4{ii},errF4{ii},timF4,xF4{ii},dxF4(ii),nF4(ii),NtotF4(ii),~] = fetchOutputs(F4);
    disp([int2str(ii),' F4 ', num2str(timF4)])
    disp([NtotF4(ii),d,p]);
    erF4(ii) = norm(errF4{ii},Inf);
    disp(erF4(ii));
    
    f=figure(runnumber);
    clf
    loglog(dxF0,erF0,'-*')
    hold on
    loglog(dxF1,erF1,'-*')
    loglog(dxF2,erF2,'-*')
    loglog(dxF3,erF3,'-*')
    loglog(dxF4,erF4,'-*')
    title('BSeuCall2Dbasket')
    legend('FD2\_reg\_square','PHSd5p5\_reg\_triangle', 'PHSd7p7\_reg\_triangle','PHSd5p5\_reg\_triangle\_smooth', 'PHSd7p7\_reg\_triangle\_smooth')
    xlabel('h')
    ylabel('\Deltau_{max}')
    drawnow
end

% disp([N,d,p]);
% disp(norm(err,Inf));
% disp(tim);

% xvec = x(:,1); yvec = x(:,2);

% figure(1)
% clf
% plot(xvec,yvec,'.')
% hold on
% plot(xvec(indff),yvec(indff),'*')
% plot(xvec(indcf),yvec(indcf),'^')
% plot(xvec(indin),yvec(indin),'o')
% axis equal
% axis tight
% hold off

% figure(1)
% tri = delaunay(xvec,yvec);
% trisurf(tri, xvec, yvec, u);
% shading interp
% colorbar
% view(2)
% axis vis3d
% axis tight
% 
% figure(2)
% % tri = delaunay(xvec,yvec);
% trisurf(tri, xvec, yvec, err);
% shading interp
% colorbar
% view(2)
% axis vis3d
% axis tight

% N = [25, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600];
%
% Kx = 1;

% %% Run
% for ii = 1:numel(N)
%
%     n=5;
%     F0 = parfeval(@BSeuCall1D_RBFFDadapsmooth,7,N(ii),n,M,Kmul);
%
%     F1 = parfeval(@BSeuCall1D_FD2half,8,N(ii),M,Kmul);
%     F2 = parfeval(@BSeuCall1D_FD4half,8,N(ii),M,Kmul);
%     F3 = parfeval(@BSeuCall1D_FD2,7,N(ii),M,Kmul);
%
%
%
%
%     d = 5;
%     p = 5;
% %     F4 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
%     F4 = parfeval(@BSeuCall1D_RBFFDreg_phs_smooth,7,N(ii),p,d,M,Kmul);
%
% %     F5 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
%     F5 = parfeval(@BSeuCall1D_RBFFDreg_phs_adap_smooth,7,N(ii),p,d,M,Kmul);
%
%     p = 7;
% %     F6 = parfeval(@BSeuCall1D_RBFFDreg_phs,7,N(ii),n,r,M,Kmul);
%     F6 = parfeval(@BSeuCall1D_RBFFDreg_phs_adap_smooth,7,N(ii),p,d,M,Kmul);
%
%     [uF0{ii},errF0{ii},timF0,xF0{ii},dxF0(ii),~,~] = fetchOutputs(F0);
%     disp([int2str(ii),' F0 ', num2str(timF0)])
%     [uF1{ii},errF1{ii},timF1,xF1{ii},dxF1(ii),~,~,~] = fetchOutputs(F1);
%     disp([int2str(ii),' F1 ', num2str(timF1)])
%     [uF2{ii},errF2{ii},timF2,xF2{ii},dxF2(ii),~,~,~] = fetchOutputs(F2);
%     disp([int2str(ii),' F2 ', num2str(timF2)])
%     [uF3{ii},errF3{ii},timF3,xF3{ii},dxF3(ii),~,~] = fetchOutputs(F3);
%     disp([int2str(ii),' F3 ', num2str(timF3)])
%
%     [uF4{ii},errF4{ii},timF4,xF4{ii},dxF4(ii),~,~] = fetchOutputs(F4);
%     disp([int2str(ii),' F4 ', num2str(timF4)])
%
%     [uF5{ii},errF5{ii},timF5,xF5{ii},dxF5(ii),~,~] = fetchOutputs(F5);
%     disp([int2str(ii),' F5 ', num2str(timF5)])
%
%     [uF6{ii},errF6{ii},timF6,xF6{ii},dxF6(ii),~,~] = fetchOutputs(F6);
%     disp([int2str(ii),' F6 ', num2str(timF6)])
%
%
%     indreg=[];
%     for jj = 1:length(xF1{ii})
%         %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%         if xF1{ii}(jj) >= 1/3*Kx && xF1{ii}(jj) <= 5/3*Kx
%             indreg=[indreg jj];
%         end
%     end
%
%     xind = xF1{ii}(indreg);
%
%     errF5crop{ii} = errF5{ii}(indreg);
%     errF6crop{ii} = errF6{ii}(indreg);
%     errF0crop{ii} = errF0{ii}(indreg);
%     errF4crop{ii} = errF4{ii}(indreg);
%     errF1crop{ii} = errF1{ii}(indreg);
%     errF2crop{ii} = errF2{ii}(indreg);
%     errF3crop{ii} = errF3{ii}(indreg);
%
%     erF5(ii) = max(abs(errF5crop{ii}));
%     erF6(ii) = max(abs(errF6crop{ii}));
%     erF0(ii) = max(abs(errF0crop{ii}));
%     erF4(ii) = max(abs(errF4crop{ii}));
%     erF1(ii) = max(abs(errF1crop{ii}));
%     erF2(ii) = max(abs(errF2crop{ii}));
%     erF3(ii) = max(abs(errF3crop{ii}));
%
%     figure(runnumber)
%     clf
%     loglog(dxF0,erF0,'-*')
%     hold on
%     loglog(dxF1,erF1,'-*')
%     loglog(dxF2,erF2,'-*')
%     loglog(dxF3,erF3,'-*')
%     loglog(dxF4,erF4,'-*')
%     loglog(dxF5,erF5,'-*')
%     loglog(dxF6,erF6,'-*')
%     legend('F0', 'F1', 'F2', 'F3', 'F4','F5', 'F6');
%     title(['M=',num2str(M)]);
%     %     xlim([1e-4, 1e-1]);
% %         ylim([1e-9, 1e-1]);
%     drawnow
% end

% %% Data
% s = ['run_M',num2str(M),'Nmax',num2str(max(N)),'K',num2str(Kmul)];
% save(s);
%
% % %% Convergence
% orderFD2 = polyfit(log10(dxF1(1:end)),log10(erF1(1:end)),1)
% orderFD4 = polyfit(log10(dxF2(1:end)),log10(erF2(1:end)),1)
% orderF3 = polyfit(log10(dxF3(1:end)),log10(erF3(1:end)),1)
% orderF0 = polyfit(log10(dxF0(1:end)),log10(erF0(1:end)),1)
% orderF4 = polyfit(log10(dxF4(1:end)),log10(erF4(1:end)),1)
% orderF5 = polyfit(log10(dxF5(1:end)),log10(erF5(1:end)),1)
% orderF6 = polyfit(log10(dxF6(1:end)),log10(erF6(1:end)),1)
%
% %% Plots
% figure()
% loglog(dxF1,erF1,'-*')
% hold on
% loglog(dxF2,erF2,'-*')
% loglog(dxF3,erF3,'-*')
% legend('F1','F2','F3ed');
%
% figure(2)
% plot(xF1{end-2},errF1{end-2})
% hold on
% plot(xF1{end-1},errF1{end-1})
% plot(xF1{end},errF1{end})
% title('T1\_F1')
%
% figure(3)
% plot(xF2{end-2},errF2{end-2})
% hold on
% plot(xF2{end-1},errF2{end-1})
% plot(xF2{end},errF2{end})
% title('T1\_F2')
%
% figure(5)
% plot(xF3{end-2},errF3{end-2})
% hold on
% plot(xF3{end-1},errF3{end-1})
% plot(xF3{end},errF3{end})
% plot(xF3{1},errF3{1})
% plot(xF3{2},errF3{2})
% plot(xF3{3},errF3{3})
% title('T1\_F3')