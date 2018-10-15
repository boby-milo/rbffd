function [W,hloc] = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,x,y,N,n,m,p,indin,phi,ep,adap,parallel,Kx)
% Constructs a BS differentiation matrix W from 2D grid s,
%stencil size n and indices indin.

if nargin == 12
    adap = 'reg';
    parallel = 0;
elseif nargin == 11
    parallel = 0;
end

% if parallel
%     ppar = gcp();
%     argfor = ppar.NumWorkers;
% else
%     ppar = gcp('nocreate');
% %     delete(p);
%     argfor = 0;
% end

indc = findKNearestNeighbors(s,s,n);

iind = repmat(indin,n,1); iind = iind(:); %n*N
%pause
jind = transpose(indc(indin,:)); jind = jind(:);%n*N
Wval = zeros(n,numel(indin));  %n*N

% parfor (ii = indin, argfor)
hloc = zeros(size(s,1),1); 
br = 0;
Ny=sqrt(N);
Nx=Ny;
 xvec=s(:,1);
 yvec=s(:,2);
%pause
dx=y(2)-y(1);
dy=dx;

% FD4

%W=spalloc(Ny*Nx,Ny*Nx,25*Ny*Nx);


%     aax = -(1/12)*(r*x)/dx +(1/24)*(sig1^2*x.^2)/dx^2;
%     bbx = (2/3)*(r*x)/dx -(2/3)*(sig1^2*x.^2)/dx^2;
%     ccx = (5/4)*(sig1^2*x.^2)/dx^2;
%     ddx = -(2/3)*(r*x)/dx -(2/3)*(sig1^2*x.^2)/dx^2;
%     eex = (1/12)*(r*x)/dx +(1/24)*(sig1^2*x.^2)/dx^2;
%     
%     aax(1)=0;
%     bbx(1)=0;
%     ccx(1)=0;
%     ddx(1)=0;
%     eex(1)=0;
%     
%     aax(2)=0;
%     bbx(2)=(1/2)*(r*dx)/dx-(1/2)*(sig1^2*(dx)^2)/dx^2;
%     ccx(2)=1*(sig1^2*(dx)^2)/dx^2;
%     ddx(2)=-(1/2)*(r*dx)/dx-(1/2)*(sig1^2*(dx)^2)/dx^2;
%     eex(2)=0;
%     
%     aax(Nx-1)=0;
%     bbx(Nx-1)=(1/2)*(r*(1-dx))/dx-(1/2)*(sig1^2*(1-dx)^2)/dx^2;
%     ccx(Nx-1)=1*(sig1^2*(1-dx)^2)/dx^2;
%     ddx(Nx-1)=-(1/2)*(r*(1-dx))/dx-(1/2)*(sig1^2*(1-dx)^2)/dx^2;
%     eex(Nx-1)=0;
%     
%     aay = -(1/12)*(r*y)/dy +(1/24)*(sig2^2*y.^2)/dy^2;
%     bby = (2/3)*(r*y)/dy -(2/3)*(sig2^2*y.^2)/dy^2;
%     ccy = (5/4)*(sig2^2*y.^2)/dy^2;
%     ddy = -(2/3)*(r*y)/dy -(2/3)*(sig2^2*y.^2)/dy^2;
%     eey = (1/12)*(r*y)/dy +(1/24)*(sig2^2*y.^2)/dy^2;
%     
%     aay(1)=0;
%     bby(1)=0;
%     ccy(1)=0;
%     ddy(1)=0;
%     eey(1)=0;
%     
%     aay(2)=0;
%     bby(2)=(1/2)*(r*dy)/dy-(1/2)*(sig2^2*(dy)^2)/dy^2;
%     ccy(2)=1*(sig2^2*(dy)^2)/dy^2;
%     ddy(2)=-(1/2)*(r*dy)/dy-(1/2)*(sig2^2*(dy)^2)/dy^2;
%     eey(2)=0;
%     
%     aay(Ny-1)=0;
%     bby(Ny-1)=(1/2)*(r*(1-dy))/dy-(1/2)*(sig2^2*(1-dy)^2)/dy^2;
%     ccy(Ny-1)=1*(sig2^2*(1-dy)^2)/dy^2;
%     ddy(Ny-1)=-(1/2)*(r*(1-dy))/dy-(1/2)*(sig2^2*(1-dy)^2)/dy^2;
%     eey(Ny-1)=0;
%     
%     aaax = -(1/12)*(rho*sig1*x)/dx;
%     bbbx = (2/3)*(rho*sig1*x)/dx;
%     cccx = zeros(size(bbbx));
%     dddx = -(2/3)*(rho*sig1*x)/dx ;
%     eeex = (1/12)*(rho*sig1*x)/dx ;
%     
%     aaax(1)=0;
%     bbbx(1)=0;
%     cccx(1)=0;
%     dddx(1)=0;
%     eeex(1)=0;
%     
%     aaax(2)=0;
%     bbbx(2)=(1/2)*(rho*sig1*dx)/dx;
%     cccx(2)=0;
%     dddx(2)=-(1/2)*(rho*sig1*dx)/dx;
%     eeex(2)=0;
%     
%     aaax(Nx-1)=0;
%     bbbx(Nx-1)=(1/2)*(rho*sig1*(1-dx))/dx;
%     cccx(Nx-1)=0;
%     dddx(Nx-1)=-(1/2)*(rho*sig1*(1-dx))/dx;
%     eeex(Nx-1)=0;
%     
%     aaay = -(1/12)*(sig2*y)/dy;
%     bbby = (2/3)*(sig2*y)/dy;
%     cccy = zeros(size(bbby));
%     dddy = -(2/3)*(sig2*y)/dy;
%     eeey = (1/12)*(sig2*y)/dy;
%     
%     aaay(1)=0;
%     bbby(1)=0;
%     cccy(1)=0;
%     dddy(1)=0;
%     eeey(1)=0;
%     
%     aaay(2)=0;
%     bbby(2)=(1/2)*(sig2*dy)/dy;
%     cccy(2)=0;
%     dddy(2)=-(1/2)*(sig2*dy)/dy;
%     eeey(2)=0;
%     
%     aaay(Ny-1)=0;
%     bbby(Ny-1)=(1/2)*(sig2*(1-dy))/dy;
%     cccy(Ny-1)=0;
%     dddy(Ny-1)=-(1/2)*(sig2*(1-dy))/dy;
%     eeey(Ny-1)=0;
%    
%     
%     j=1;
%     
% %     k=1; 
% %     
% %     W((j-1)*Ny+k,(j-1)*Ny+k)=ccx(j)+ccy(k)+r;
% %     W((j-1)*Ny+k+1,(j-1)*Ny+k)=ddx(j);
% %     W((j-1)*Ny+k+2,(j-1)*Ny+k)=eex(j);
% %     W((j-1)*Ny+Ny+k,(j-1)*Ny+k)=ddy(k);
% %     W((j-1)*Ny+2*Ny+k,(j-1)*Ny+k)=eey(k);
%     
%     k=2;
% 
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     for k=3:Ny-3
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     end
%     
%     k=Ny-2;
% 
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     k=Ny-1;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%         
%     j=2;
%     
%     k=1;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+k+1)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k+2)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     k=2;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k-1)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     for k=3:Ny-3
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     end
%     
%     k=Ny-2;
% 
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     k=Ny-1;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     for j=3:Nx-3
%         k=1;
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+k+1)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k+2)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%         k=2;
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+k-1)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%         for k=3:Ny-3
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%             W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%             W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%             W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%             W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%             W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%             W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%             
%             
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny-2+k)=-aaax(j)*aaay(k);
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny-1+k)=-aaax(j)*bbby(k);
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)-aaax(j)*cccy(k);
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny+1+k)=-aaax(j)*dddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny-2*Ny+2+k)=-aaax(j)*eeey(k);
%                     
%             W((j-1)*Ny+k,(j-1)*Ny-Ny-2+k)=-bbbx(j)*aaay(k);
%             W((j-1)*Ny+k,(j-1)*Ny-Ny-1+k)=-bbbx(j)*bbby(k);
%             W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=W((j-1)*Ny+k,(j-1)*Ny-Ny+k)-bbbx(j)*cccy(k);
%             W((j-1)*Ny+k,(j-1)*Ny-Ny+1+k)=-bbbx(j)*dddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny-Ny+2+k)=-bbbx(j)*eeey(k);
%             
%             W((j-1)*Ny+k,(j-1)*Ny-2+k)=W((j-1)*Ny+k,(j-1)*Ny-2+k)-cccx(j)*aaay(k);
%             W((j-1)*Ny+k,(j-1)*Ny-1+k)=W((j-1)*Ny+k,(j-1)*Ny-1+k)-cccx(j)*bbby(k);
%             W((j-1)*Ny+k,(j-1)*Ny+k)=W((j-1)*Ny+k,(j-1)*Ny+k)-cccx(j)*cccy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+1+k)=W((j-1)*Ny+k,(j-1)*Ny+1+k)-cccx(j)*dddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2+k)=W((j-1)*Ny+k,(j-1)*Ny+2+k)-cccx(j)*eeey(k);
%             
%             W((j-1)*Ny+k,(j-1)*Ny+Ny-2+k)=-dddx(j)*aaay(k);
%             W((j-1)*Ny+k,(j-1)*Ny+Ny-1+k)=-dddx(j)*bbby(k);
%             W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=W((j-1)*Ny+k,(j-1)*Ny+Ny+k)-dddx(j)*cccy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+Ny+1+k)=-dddx(j)*dddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+Ny+2+k)=-dddx(j)*eeey(k);
%             
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny-2+k)=-eeex(j)*aaay(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny-1+k)=-eeex(j)*bbby(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)-eeex(j)*cccy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny+1+k)=-eeex(j)*dddy(k);
%             W((j-1)*Ny+k,(j-1)*Ny+2*Ny+2+k)=-eeex(j)*eeey(k);
% 
%         end
%         k=Ny-2;
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%         k=Ny-1;
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     end
%     
%     j=Nx-2;
%     
%     k=1;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+k+1)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k+2)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     k=2;
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k-1)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     for k=3:Ny-3
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%         W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     end
%     
%     k=Ny-2;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     k=Ny-1;
%     
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+2*Ny+k)=eex(j);
%     
%     j=Nx-1;
%     
%     k=1;
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+k+1)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k+2)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
%  
%     k=2;
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny+k-1)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
% 
%     for k=3:Ny-3
%         W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%         W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%         W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%         W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%         W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%         W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%         W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%         W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
% 
%     end
%     
%     k=Ny-2;
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+2+k)=eey(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
% 
%     k=Ny-1;
%     W((j-1)*Ny+k,(j-1)*Ny-2*Ny+k)=aax(j);
%     W((j-1)*Ny+k,(j-1)*Ny-Ny+k)=bbx(j);
%     W((j-1)*Ny+k,(j-1)*Ny-2+k)=aay(k);
%     W((j-1)*Ny+k,(j-1)*Ny-1+k)=bby(k);
%     W((j-1)*Ny+k,(j-1)*Ny+k)=ccy(k)+ccx(j)+r;
%     W((j-1)*Ny+k,(j-1)*Ny+1+k)=ddy(k);
%     W((j-1)*Ny+k,(j-1)*Ny+Ny+k)=ddx(j);
% 
%     W=-W;
%  
% RBF-FD

    for ii=[indin]
    br = br + 1;
    sc = s(ii,:); xc = sc(:,1); yc = sc(:,2);
    se = s(indc(ii,:),:);
    
    Rc = xcdist(se,se,1);
    H = Rc(:,:,1);
    hloc(ii) = min(H(H>0));
    
%     if strcmp(adap,'min')
%         H = Rc(:,:,1);
%         hmin = min(H(H>0));
% %         ep = gamma/hmin;
%         
%     elseif strcmp(adap,'avg')
%         H = Rc(:,:,1);
%         havg = mean(H(H>0));
% %         ep = gamma/havg;
%     end
    
    A = RBFmat(phi,ep,Rc,'0',1);
    
    Ax = RBFmat(phi,ep,Rc,'1',1);
    Ay = RBFmat(phi,ep,Rc,'1',2);
    
    Axx = RBFmat(phi,ep,Rc,'2',1);
    Ayy = RBFmat(phi,ep,Rc,'2',2);
    Axy = RBFmat(phi,ep,Rc,'m2',1:2);
    
    [P,lP] = vander2D(sc,se,n,p,r,sig1,sig2,rho);
    
    
    Ac = [A, P;
          P', zeros(m,m)];
    
    lA = transpose(-r*A(1,:)...
        +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)...
        +0.5*sig1^2*xc'.^2.*Axx(1,:)...
        +0.5*sig2^2*yc'.^2.*Ayy(1,:)...
        +rho*sig1*sig2*xc'.*yc'.*Axy(1,:));

    
    
   lc = [lA; lP];
    
    wc = Ac\lc;
    Wval(:,br) = wc(1:n);
    end
%         end
%     end
% end

Wval = Wval(:);
W = sparse(iind,jind,Wval,N,N);%-Wp;

% for ii=1:N
%     if (sum(ismember(ii,mynoll))>0)
% %         ii
% %         pause
%         W(:,ii)=zeros(N,1);
%     end
% end
% figure(222)
% spy(W)
% figure(333)
% spy(Wp)
% 
% W=W-Wp;
% figure(444)
% spy(W)
% pause
%     end
% end
end


function [P, lP] = vander2D(sc,s,n,p,r,sig1,sig2,rho)
x = s(:,1);
y = s(:,2);

xc = sc(1);
yc = sc(2);

dim = 2;
m = nchoosek(p+dim, p); %number of polynomial terms;

P = zeros(n,m);
lP = zeros(m,1);

kk = 0;
for ii = 0:p
    for jj = 0:ii
        kk = kk+1;
        
        i = ii-jj;
        j = jj;
        
        P(:,kk) = x.^(i) .* y.^(j);
        
        lP(kk) = xc^(i)*yc^(j) * (-r +r*(i+j) + 0.5*(i-1)*i*sig1^2 + 0.5*(j-1)*j*sig2^2 +i*j*rho*sig1*sig2);
        
    end
end
end