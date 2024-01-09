function [rel_err,err]=normQ3d(T,Tap)
sTap3=tens2mat(s(Tap),3);
xTap3=tens2mat(x(Tap),3);
yTap3=tens2mat(y(Tap),3);
zTap3=tens2mat(z(Tap),3);

sT3=tens2mat(s(T),3);
xT3=tens2mat(x(T),3);
yT3=tens2mat(y(T),3);
zT3=tens2mat(z(T),3);

MaxTap=quaternion(sTap3,xTap3,yTap3,zTap3);
MaxT=quaternion(sT3,xT3,yT3,zT3);

 
difT0_compress=MaxT-MaxTap;
err=norm(difT0_compress,'fro');
rel_err=err/norm(MaxT,'fro');