function [phi]=phs(m,r,nprime,dim)
% This is a polyharmonic spline (PHS) phi(r) = r^m.
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% DIM(1:2) are the dimensions for the mixed second derivative if nprime
% is 'm2'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
%
% Assume a one-dimensional problem if no dimension is given
%
if nargin <= 3
    dim = 1;
end
%
% For the mixed derivative, the dimensions must be given even in 2D
%
if (nprime(1)=='m')
    if (length(dim)~=2)
        error('For the mixed derivative, dim=dim(1:2)')
    elseif (dim(1)==dim(2))
        error('For mixed derivatives, dim(1) must be other than dim(2)')
    end
end

if nprime(1)=='0'
    phi = sq(r(:,:,1)).^m;
    
elseif nprime(1)=='1'
    phi = m*sq(r(:,:,1).^(m-2).*r(:,:,dim+1));
    
elseif nprime(1)=='2'
    % make sure that phi is computed correctly for r=0
    mask=(r(:,:,1)==0);
    phi = m*(m-2)*sq((r(:,:,1)+mask).^(m-4).*r(:,:,dim+1).^2) + m*sq(r(:,:,1).^(m-2));
    
elseif nprime(1:2)=='m2'
    % make sure that phi is computed correctly for r=0
    mask=(r(:,:,1)==0);
    phi = m*(m-2)*sq((r(:,:,1)+mask).^(m-4).*r(:,:,dim(1)+1).*r(:,:,dim(2)+1));
    
    % elseif nprime(1)=='L' & length(nprime)==1
    % nd = size(r,3)-1;
    % phi = 3*(nd+1)*sq(r(:,:,1));
    
else
    error('Error in input argument nprime to function phs')
end

function r=sq(r)
r=squeeze(r);
