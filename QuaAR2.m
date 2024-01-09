function [X,y0]=QuaAR2(x1,x2,toa,tob,eta,m)
% This function generate sequence by AR2 model: X(t)=-toa*X(t-1)-tob*X(t-2?+v(t),
% where v(t) is a white noise process with variance eta^2.
% m = 500 ; % sample number
% eta=0.025; toa=0.9; tob=0.5; x1=randq(1);
% x1=randq(1); x2=randq(1); 

%% ======================== generate model by AR (2) ====================

addpath('qtfm')

%===== generate X0=[x1,x2,...xm,x_{m+1}]======= 

X0=zerosq(m+1,1); X0(1)=x1; X0(2)=x2;
 
for t=3:m+1
svt=eta*randn(1); ivt=eta*randn(1); jvt=eta*randn(1); kvt=eta*randn(1);
vt=quaternion(svt,ivt,jvt,kvt);
X0(t)=-toa*X0(t-1)-tob*X0(t-2)+vt;
end 

X=X0(1:m);  y0=X0(end:-1:2);
