clear all
clc
close all

path='./1Da2g5';

%% FD
filepaths=getfilenames(path,'BSeuCall1D_FD*.mat');
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

figure(1)
loglog(hvec,errvec,'k-*');
hold on

figure(2)
loglog(timvec,errvec,'k-*');
hold on

%% RBF-FD const
% n=3
filepaths=getfilenames(path,'BSeuCall1D_RBFFDc1reg_ep3*.mat');
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
figure(1)
loglog(hvec,errvec,'r*:');

figure(2)
loglog(timvec,errvec,'r*:');

% n=5
filepaths=getfilenames(path,'BSeuCall1D_RBFFDc1reg_ep5*.mat');
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
figure(1)
loglog(hvec,errvec,'b*:');


figure(2)
loglog(timvec,errvec,'b*:');



% n=7
filepaths=getfilenames(path,'BSeuCall1D_RBFFDc1reg_ep7*.mat');
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


figure(1)
loglog(hvec,errvec,'g*:');
% ylim([1e-6 1e-1])
xlim([1e-2 1e-0])
axis tight

figure(2)
loglog(timvec,errvec,'g*:');

%% RBF-FD adap
% n=3
filepaths=getfilenames(path,'BSeuCall1D_RBFFDadapEps_ep3*.mat');
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
figure(1)
loglog(hvec,errvec,'r*-');

figure(2)
loglog(timvec,errvec,'r*-');

% n=5
filepaths=getfilenames(path,'BSeuCall1D_RBFFDadapEps_ep5*.mat');
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
figure(1)
loglog(hvec,errvec,'b*-');


figure(2)
loglog(timvec,errvec,'b*-');



% n=7
filepaths=getfilenames(path,'BSeuCall1D_RBFFDadapEps_ep7*.mat');
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


figure(1)
loglog(hvec,errvec,'g*-');
% ylim([1e-6 1e-1])
xlim([1e-2 1e-0])
axis tight

figure(2)
loglog(timvec,errvec,'g*-');

%% RBF-FD reg
% n=3
filepaths=getfilenames(path,'BSeuCall1D_RBFFDreg_ep3*.mat');
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
figure(1)
loglog(hvec,errvec,'r-.*');

figure(2)
loglog(timvec,errvec,'r-.*');

% n=5
filepaths=getfilenames(path,'BSeuCall1D_RBFFDreg_ep5*.mat');
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
figure(1)
loglog(hvec,errvec,'b-.*');


figure(2)
loglog(timvec,errvec,'b-.*');



% n=7
filepaths=getfilenames(path,'BSeuCall1D_RBFFDreg_ep7*.mat');
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


figure(1)
loglog(hvec,errvec,'g-.*');
ylim([1e-6 1e-1])
xlim([1e-2 1e-0])
axis tight
legend('FD','con3','con5','con7','adap3','adap5','adap7','reg3','reg5','reg7')
title('European Call 1D')
xlabel('h')
ylabel('\Deltau')
grid on

figure(2)
loglog(timvec,errvec,'g-.*');
ylim([1e-6 1e-1])
% xlim([1e-2 0.3])
axis tight
legend('FD','con3','con5','con7','adap3','adap5','adap7','reg3','reg5','reg7')
title('European Call 1D')
xlabel('t')
ylabel('\Deltau')
grid on


