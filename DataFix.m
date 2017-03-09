clc
close all
clear

%% FD

prename = 'phs_test_';
basename = 'BSeuCall2Dbasket_FD__';

fullname = [prename, basename];

funname = 'BSeuCall2Dbasket_FD';


names = dir(['./Data/',fullname,'*']);

for ii = 1:numel(names)
    cd('./Data')
    load(names(ii).name)
    N = Ntot;
    Nx = Nin;
    M = Min;
    
    datafilename = [funname,'___Nx',num2str(Nx),'_M',num2str(M),'_Kmul',num2str(Kmul),'2.mat'];
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'N', 'W', 'Nx', 'M', 'Kmul')
    cd('..')
end

%% RBFFDadap
prename = 'phs_test_';
basename = 'BSeuCall2Dbasket_RBFFDadap_phs___';

fullname = [prename, basename];

funname = 'BSeuCall2Dbasket_RBFFDadap_phs';


names = dir(['./Data/',fullname,'*']);

for ii = 1:numel(names)
    cd('./Data')
    load(names(ii).name)
    N = Ntot;
    Nx = Nin;
    M = Min;
    
    datafilename = [funname,'___Nx',num2str(Nx),'_g',num2str(g),'_p',num2str(p),'_d',num2str(d),'_nm',num2str(nm),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'N', 'W', 'Nx', 'M', 'Kmul')
    cd('..')
end
