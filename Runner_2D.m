clc
clear all
close all

warning off

cd('./2Da2fit1.75new')
a=2;
fit=2;

gamma9=0.011; gamma9=fit*gamma9;
gamma13=0.0145; gamma13=fit*gamma13;
gamma25=0.045; gamma25=fit*gamma25;

% h=linspace(0.05,0.8,100];
i=linspace(0.5,2,10);
Nx=round(8*10.^i)
hx=8./(Nx-1)


M=round(a./hx)

eps9=gamma9./hx; eps9c1=1.3; eps9c2=0.05;
eps13=gamma13./hx; eps13c1=4; eps13c2=0.1;
eps25=gamma25./hx; eps25c1=3; eps25c2=0.3;



for ii=1:numel(Nx)
    display([num2str(ii/numel(Nx)*100),'%'],'Progress')
    
    %% FD
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % EU
    [u,err,tim,s,h,N,~]=BSeuCallbasket2D_FD(Nx(ii),M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time FD EUCallbasket')
    filename=['BSeuCallbasket2D_FD_','Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % AM
    [u,err,tim,s,h,N,~]=BSamPutbasket2D_FD(Nx(ii),M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time FD AMPutbasket')
    filename=['BSamPutbasket2D_FD_','Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M');
    
    
    
    %% RBF-FD n=9
    eps=eps9(ii)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % EU   
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDreg(Nx(ii),9,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg EUCallbasket')
    filename=['BSeuCallbasket2D_RBFFDreg_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsmin(Nx(ii),9,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time EU adapEpsmin')
    filename=['BSeuCallbasket2D_RBFFDadapEpsmin_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsavg(Nx(ii),9,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time EU adapEpsavg')
    filename=['BSeuCallbasket2D_RBFFDadapEpsavg_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % AM
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDreg(Nx(ii),9,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg AMPutbasket')
    filename=['BSamPutbasket2D_RBFFDreg_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsmin(Nx(ii),9,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time AM adapEpsmin')
    filename=['BSamPutbasket2D_RBFFDadapEpsmin_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsavg(Nx(ii),9,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time AM adapEpsavg')
    filename=['BSamPutbasket2D_RBFFDadapEpsavg_','ep',num2str(9),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    
    
    
    %% RBF-FD n=13
    eps=eps13(ii)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % EU
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDreg(Nx(ii),13,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg EUCallbasket')
    filename=['BSeuCallbasket2D_RBFFDreg_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsmin(Nx(ii),13,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time EU adapEpsmin')
    filename=['BSeuCallbasket2D_RBFFDadapEpsmin_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
     [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsavg(Nx(ii),13,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time EU adapEpsavg')
    filename=['BSeuCallbasket2D_RBFFDadapEpsavg_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % AM
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDreg(Nx(ii),13,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg AMPutbasket')
    filename=['BSamPutbasket2D_RBFFDreg_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsmin(Nx(ii),13,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time AM adapEpsmin')
    filename=['BSamPutbasket2D_RBFFDadapEpsmin_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsavg(Nx(ii),13,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time AM adapEpsavg')
    filename=['BSamPutbasket2D_RBFFDadapEpsavg_','ep',num2str(13),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    
     %% RBF-FD n=25
    eps=eps25(ii)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % EU
    
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDreg(Nx(ii),25,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg EUCallbasket')
    filename=['BSeuCallbasket2D_RBFFDreg_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsmin(Nx(ii),25,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time EU adapEpsmin')
    filename=['BSeuCallbasket2D_RBFFDadapEpsmin_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
     [u,err,tim,s,h,N,~] = BSeuCallbasket2D_RBFFDadapEpsavg(Nx(ii),25,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time EU adapEpsavg')
    filename=['BSeuCallbasket2D_RBFFDadapEpsavg_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % AM
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDreg(Nx(ii),25,eps,M(ii));
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h),', E = ', num2str(max(abs(err)))],'Time RBFFDreg AMPutbasket')
    filename=['BSamPutbasket2D_RBFFDreg_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsmin(Nx(ii),25,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time AM adapEpsmin')
    filename=['BSamPutbasket2D_RBFFDadapEpsmin_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
    [u,err,tim,s,h,N,~] = BSamPutbasket2D_RBFFDadapEpsavg(Nx(ii),25,M(ii),fit);
    display([num2str(tim),' s for Nx = ',num2str(Nx(ii)),', h = ',num2str(h), ', E = ', num2str(max(abs(err)))],'Time AM adapEpsavg')
    filename=['BSamPutbasket2D_RBFFDadapEpsavg_','ep',num2str(25),'Nx',num2str(Nx(ii))];
    save(filename,'u','err','tim','s','h','N','M','eps');
    
end





