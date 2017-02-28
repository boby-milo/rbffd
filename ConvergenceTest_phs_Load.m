ConvergenceTest_phs

%% ConvergenceTest_phs

%FD
names = dir('./Data/phs*FD__*');
names.name
legplot{1} = 'FD2';
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
    
end
[h,hsort] = sort(h);
delu = delu(hsort);

xo = 10.^(-2:0.1:0);
yo2 = (1/10)*xo.^2;
yo4 = (1/10)*xo.^4;
figure(1)
loglog(h,delu,'-*')
hold on

%RBFFDreg
names = dir('./Data/phs*RBFFDreg_phs*N*p7*d7*nm3*');
names.name
legplot{2} = 'RBFFDregp7d7nm3';
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
end
[h,hsort] = sort(h);
delu = delu(hsort);
figure(1)
loglog(h,delu,'-o')

% RBFFDadap
names = dir('./Data/phs*RBFFDadap_phs*N*p5*d5*nm3*g5*');
names.name
legplot{3} = 'RBFFDadapp5d5nm3g5';
for ii = 1:numel(names)
    
    cd('./Data')
    load(names(ii).name);
    cd('..')
    
    delu(ii) = max(abs(err));
    h(ii) = dx;
end
[h,hsort] = sort(h);
delu = delu(hsort);
figure(1)
loglog(h,delu,'-^')
legend(legplot)
loglog(xo,yo2,'k--')
loglog(xo,yo4,'k--')
xlim([1e-2,1e0]);
ylim([1e-9,1e-1]);

% names = dir('./Data/phs*RBFFDreg_smooth*N40*nm1.5*')
% for ii = 1:numel(names)
%     
%     cd('./Data')
%     load(names(ii).name);
%     cd('..')
%     
%     delu(ii) = max(abs(err));
%     h(ii) = dx;
%     legname{ii} = names(ii).name;
% end
% delu'
% legname'

