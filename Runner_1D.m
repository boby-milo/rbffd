clc
clear all
close all

warning off

cd('./1Da2g5fit1.2')

% 1D
gamma3=10^(-3); gamma3=1.2*gamma3;
gamma5=0.0145; gamma5=1.2*gamma5;
gamma7=0.045; gamma7=1.2*gamma7;

i=linspace(1,4,10);
Nx=round(4*10.^i)
hx=4./(Nx-1)

a=2;
M=round(a./hx);

eps3=gamma3./hx; eps3c1=3.02; eps3c2=0.01;
eps5=gamma5./hx; eps5c1=43.81; eps5c2=0.5;
eps7=gamma7./hx; eps7c1=135.98; eps7c2=1.5;

fit=1.2;

for ii=1:numel(Nx)
    display([num2str(ii/numel(Nx)*100),'%'],'Progress')
    
    %% FD
    
    [u,err,tim,s,h,N,~]=BSeuCall1D_FD(Nx(ii),M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time FD')
    filename=['BSeuCall1D_FD_','Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M');    
    
    %% RBF-FD n=3
    eps=eps3(ii)
    
    %epsconst
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),3,eps3c1,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDc1reg_','ep',num2str(3),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),3,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDreg_','ep',num2str(3),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
        
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDadapEps(Nx(ii),3,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDadap')
    filename=['BSeuCall1D_RBFFDadapEps_','ep',num2str(3),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');    
    
    %% RBF-FD n=5
    eps=eps5(ii)

    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),5,eps5c1,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDc1reg_','ep',num2str(5),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),5,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDreg_','ep',num2str(5),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
        
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDadapEps(Nx(ii),5,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDadap')
    filename=['BSeuCall1D_RBFFDadapEps_','ep',num2str(5),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');  

    
     %% RBF-FD n=7
    eps=eps7(ii)
    
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),7,eps7c1,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDc1reg_','ep',num2str(7),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDreg(Nx(ii),7,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg')
    filename=['BSeuCall1D_RBFFDreg_','ep',num2str(7),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
        
    [u,err,tim,s,h,N,~] = BSeuCall1D_RBFFDadapEps(Nx(ii),7,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDadap')
    filename=['BSeuCall1D_RBFFDadapEps_','ep',num2str(7),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');  
    
end





