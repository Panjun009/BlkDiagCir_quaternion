function c=TCpreconder(Lrvec)
% This function is generate Tony Chan's precondintoner for Hermitian
% Toeplitz matrix. 
% Lrvec=[rvec,rvec_inv];
% rvec=[r0,r1,...r_{N-1}]; rvec_inv=[r_{-1},r_{-2},...r_{-N+1},r_{-N}];
% for Hermitian Toeplitz matrix: rvec is the 1st column of Hermitian
% Toeplitz; rvec_inv is  (2nd-Nth elements of the 1st row, r(-N)) 
 
m=length(Lrvec);n=m/2;
 
rvec=Lrvec(1:n);  rvec_inv=Lrvec(n+1:m);
 
c=zerosq(n,1);
 
for t=1:n
    j=t-1;
    up=(n-j)*rvec(t)+j*rvec_inv(n-t+1);
    c(t)= up/n;     
end