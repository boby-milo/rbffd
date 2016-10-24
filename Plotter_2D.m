clear all
clc
close all

path='./2da2fit3';

a=5;
b=6;

%% FD
filepaths=getfilenames(path,'BSeuCallbasket2D_FD*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);

figure(a)
clf
loglog(hvec,errvec,'k-o');
hold on

figure(b)
clf
loglog(timvec,errvec,'k-o');
hold on

%% RBF-FD adap
% n=9
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDadapEpsmin_ep9*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);
figure(a)
loglog(hvec,errvec,'r-*');

figure(b)
loglog(timvec,errvec,'r-*');

% n=13
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDadapEpsmin_ep13*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);
figure(a)
loglog(hvec,errvec,'b-*');


figure(b)
loglog(timvec,errvec,'b-*');



% n=25
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDadapEpsmin_ep25*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);


figure(a)
loglog(hvec,errvec,'g-*');

figure(b)
loglog(timvec,errvec,'g-*');

%% RBF-FD reg
% n=9
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDreg_ep9*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);
figure(a)
loglog(hvec,errvec,'r-.*');

figure(b)
loglog(timvec,errvec,'r-.*');

% n=13
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDreg_ep13*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);
figure(a)
loglog(hvec,errvec,'b-.*');


figure(b)
loglog(timvec,errvec,'b-.*');



% n=25
filepaths=getfilenames(path,'BSeuCallbasket2D_RBFFDreg_ep25*.mat');
N=length(filepaths);
for ii=1:N
    load(filepaths{ii})
    
    hvec(ii)=h;
    errvec(ii)=max(abs(err));
    timvec(ii)=tim;
    
end
[hvec,ivec]=sort(hvec);
errvec=errvec(ivec);
timvec=timvec(ivec);


figure(a)
loglog(hvec,errvec,'g-.*');

% xlim([1e-2 1e-0])
axis tight
ylim([1e-5 1e-1])
legend('FD','adap9','adap13','adap25','reg9','reg13','reg25')
title('BSeuCallbasket2D')
xlabel('h')
ylabel('\Deltau')
grid on

figure(b)
loglog(timvec,errvec,'g-.*');
axis tight
ylim([1e-5 1e-1])
% xlim([1e-2 0.3])

legend('FD','adap9','adap13','adap25','reg9','reg13','reg25')
title('BSeuCallbasket2D')
xlabel('t')
ylabel('\Deltau')
grid on