clc
clear
dbstop if error
ConvergenceTest_phs

K = 1;

aplot = 1.5;
bplot = 8;

%% Plot
xo = 10.^(-3:0.1:0);
yo2 = 1/10*xo.^2;
yo4 = xo.^4;

fign = 12;

figure(fign)
clf
loglog(xo,yo2,'k--')
hold on
loglog(xo,yo4,'k--')
xlim([8e-3,2e-1]);
ylim([1e-7,1e-3]);

br = 2;
datafilename{1} = 'O2';
datafilename{2} = 'O4';

ploop = p;
dloop = d;
nmloop = nm;


%% FD2
datafilename{br} = ['./Data/*FD2*'];
names = dir(datafilename{br});

delu = [];
h = [];
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
    
    
end
[h,hsort] = sort(h);
delu = delu(hsort);
figure(fign)
loglog(h,delu,'-kd','LineWidth',aplot,'MarkerSize',bplot)
hold on

% legend('FD2');
drawnow
names.name

%% GS reg
datafilename{br} = ['./Data/*RBFFDreg*n25*'];
names = dir(datafilename{br});

delu = [];
h = [];
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
    
    
end
[h,hsort] = sort(h);
delu = delu(hsort);
figure(fign)
loglog(h,delu,'-sq','LineWidth',aplot,'MarkerSize',bplot)
hold on

% legend('GS reg');
drawnow
names.name

%% GS adap
datafilename{br} = ['./Data/*RBFFDadap*g5*n25*'];
names = dir(datafilename{br});

delu = [];
h = [];
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
    
    
end
[h,hsort] = sort(h);
delu = delu(hsort);
figure(fign)
loglog(h,delu,'-^','LineWidth',aplot,'MarkerSize',bplot)
hold on

%% PHS reg

for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDreg_phs_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            delu = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-sq','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
            legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end

%% PHS adap
for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDadap_phs_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            delu = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-^','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
%             legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end

%% PHS repel
for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDrepel_phs_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            del = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-o','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
%             legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end


%% PHS reg smoothed
ploop = p;
dloop = d;
nmloop = nm;

for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDreg_smoothed_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            delu = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-.sq','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
            legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end

%% PHS adap smoothed
for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDadap_smoothed_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            delu = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-.^','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
%             legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end

%% PHS repel smoothed
for ll = 1:numel(ploop)
    for jj = 1:numel(dloop)
        for kk = 1:numel(nmloop)
            br = br + 1;
            datafilename{br} = ['./Data/*RBFFDrepel_smoothed_*p',num2str(ploop(ll)),'*d',num2str(dloop(jj)),'*nm',num2str(nmloop(kk)),'*'];
            names = dir(datafilename{br});
            
            del = [];
            h = [];
            for ii = 1:numel(names)
                
                cd('./Data')
                load(names(ii).name);
                cd('..')
                
                delu(ii) = max(abs(err));
                h(ii) = dx;
                
                
            end
            [h,hsort] = sort(h);
            delu = delu(hsort);
            figure(fign)
            loglog(h,delu,'-.o','LineWidth',aplot,'MarkerSize',bplot)
            hold on
            
%             legend(datafilename);
            drawnow
            names.name
            %             pause()
            
        end
    end
end

legend('O2','O4','fd2','gs reg','gs adap','phs reg','phs adap','phs repel','phs reg smoothed','phs adap smoothed','phs repel smoothed','Location','NW')
% legend('O2','O4','gs reg','gs adap','phs reg','phs adap','phs reg smoothed','phs adap smoothed','Location','NW')
% legend('O2','O4','FD2','GS reg n25','GS adap n25', 'PHS reg p7 d7','PHS adap p7 d7','PHS reg p7 d7 smoothed','PHS adap p7 d7 smoothed','Location','NW');
% legend('O2','O4','PHS reg p7 d7','PHS adap p7 d7','PHS reg p7 d7 smoothed','PHS adap p7 d7 smoothed','Location','NW');

drawnow
names.name

xlabel('h');
ylabel('\Deltau_{max}');
title('p4d7');
set(gca,'FontSize',30)
grid on