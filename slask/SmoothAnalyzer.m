clc
clear
close all

M=1000000;
Nmax=25600;
Nmin=400;

s=['run_M',num2str(M),'Nmax',num2str(Nmax),'Nmin',num2str(Nmin),'.mat'];

load(s)

%% Convergence
orderFD2 = polyfit(log10(dxFD2half(1:4)),log10(erFD2half(1:4)),1)
orderFD4 = polyfit(log10(dxFD4half(1:3)),log10(erFD4half(1:3)),1)
orderFD4smooth = polyfit(log10(dxFD4smooth(5:6)),log10(erFD4smooth(5:6)),1)

%% Plots
figure()
loglog(dxFD2half,erFD2half,'-*')
hold on
loglog(dxFD4half,erFD4half,'-*')
loglog(dxFD4smooth,erFD4smooth,'-*')
legend('FD2half','FD4half','FD4smoothed');

% figure(2)
% plot(xFD2half{end-2},errFD2half{end-2})
% hold on
% plot(xFD2half{end-1},errFD2half{end-1})
% plot(xFD2half{end},errFD2half{end})
% title('T1\_FD2half')
% 
% figure(3)
% plot(xFD4half{end-2},errFD4half{end-2})
% hold on
% plot(xFD4half{end-1},errFD4half{end-1})
% plot(xFD4half{end},errFD4half{end})
% title('T1\_FD4half')
% 
% figure(5)
% plot(xFD4smooth{end-2},errFD4smooth{end-2})
% hold on
% plot(xFD4smooth{end-1},errFD4smooth{end-1})
% plot(xFD4smooth{end},errFD4smooth{end})
% plot(xFD4smooth{1},errFD4smooth{1})
% plot(xFD4smooth{2},errFD4smooth{2})
% plot(xFD4smooth{3},errFD4smooth{3})
% title('T1\_FD4smooth')