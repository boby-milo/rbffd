function [W] = differentiation_matrix(m,p,s)
%Discretizes spatial differential operator.
%Inputs:
%   m.type: 'FD', 'RBF-FD';
%   m.order: 2, 4;
%Outputs:
%   W: differentiation matrix;

W = [];
switch m.type
    case 'FD'
        if s.type ~= 'uniform'
            disp('FD is implemented only for uniform grids!')
            return;
        else
            switch m.order
                case 2
                    h = s.x(2)-s.x(1);
                    
                    aa = 0.5*p.sig^2*s.x.^2/h^2 -0.5*p.r*s.x/h;
                    bb = -p.sig^2*s.x.^2/h^2 -p.r;
                    cc = 0.5*p.sig^2*s.x.^2/h^2 +0.5*p.r*s.x/h;
                    
                    W = gallery('tridiag',aa(2:end),bb,cc(1:end-1));
                    W(1,:)=zeros(1,s.N);
                    W(end,:)=zeros(1,s.N);
                case 4
                    W = [];
            end
        end
    case 'RBF-FD'
        W = [];
end
end

