function [X,y0]=QuaAR1(x1,pho,eta,m)
% This function generate sequence by AR1 model: X(t)=pho*X(t-1)+v(t),
% where v(t) is a white noise process with variance eta^2.
% m = 500 ; % sample number
% eta=1; pho=0.99; x1=randq(1);
%% ======================== generate model by AR (1) ====================
addpath('qtfm')

%%===== generate X0=[x1,x2,...xm,x_{m+1}]=======
 
X0=zerosq(m+1,1); X0(1)=x1;
   
for t=2:m+1
svt=eta*randn(1); ivt=eta*randn(1); jvt=eta*randn(1); kvt=eta*randn(1);
vt=quaternion(svt,ivt,jvt,kvt);
X0(t)=pho*X0(t-1)+vt;
end 

X=X0(1:m);  y0=X0(end:-1:2);
