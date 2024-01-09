function Ch=fastQmuliply_circvec(col,h)
%% THis function computes the product of a quaternion circulant matrix (col is its first column) and  a quaternin vector (h) in an efficient way. 
%% Note Here We only consider the situation of (i,j,k). For other situation (mu, alpha,beta), one needs to transform the circulant and vector to the corresponding axis.

%% Reference: J. Pan and M. Ng, "Block Diagonalization of Quaternion Circulant Matrices with Applications "   https://arxiv.org/abs/2302.04086  

mu=quaternion(1,0,0);alpha=quaternion(0,1,0);beta=quaternion(0,0,1);
m=length(col);

zov=zeros(m,1);

tic
sc=s(col); xc=x(col); yc=y(col);zc=z(col);

% dsc=E*sc*sqrt(m);show(dsc)
% Dsc=E*circulant(sc)*E';show(Dsc)


qsc=quaternion(sc,zov,zov,zov); 
dsc0=qfft(qsc,mu,'L') ; %show(dsc0)

qxc=quaternion(xc,zov,zov,zov); 
dxc0=qfft(qxc,mu,'L') ;
qyc=quaternion(yc,zov,zov,zov); 
dyc0=qfft(qyc,mu,'L') ;
qzc=quaternion(zc,zov,zov,zov); 
dzc0=qfft(qzc,mu,'L') ;


Eh0=qfft(h,mu,'L'); % Eh=E*h*sqrt(m); % show(Eh-Eh0)

Termult1=dsc0+dxc0*mu;
Termult2=dyc0*alpha+dzc0*beta;

fEh0=[Eh0(1);flip(Eh0(2:m))];
Tersum=Termult1.*Eh0+Termult2.*fEh0;
Ch=iqfft(Tersum,mu,'L');