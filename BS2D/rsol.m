% Copyright 2014, Slobodan Milovanovic 

function u = rsol(sig, r, K, T, s)

d1=(log(s./K)+(r+0.5*sig^2)*T)/(sig*sqrt(T));
d2=d1-sig*sqrt(T);

u=0.5*s.*(1+erf(d1./sqrt(2)))-exp(-r*T)*K*0.5.*(1+erf(d2./sqrt(2)));
end