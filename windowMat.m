function T = windowMat(X,n,option)
%% X is the input samples [x1,x2,... xm], n is the window width.

m=length(X);
xL=1;

switch option
    
case 1
%% +++++++++++ by correlation method (T1) +++++++++++

d=(n-1)*xL;
Zo=zerosq(d,1); XZo=[X; Zo]; % length: xL*(m+n-1)
 
sizT=(m+n-1)*xL;
T = zerosq(sizT,n);
Xzot=zerosq(sizT,1);

for t=1:n
    if t==1
        T(:,1)= XZo;
    elseif t==n
         Xzot(1:(n-1)*xL)= zerosq(d,1);
         Xzot((n:n-1+m)*xL)=X;
         T(:,t)= Xzot;  
    else 
    Xzot( (1:t-1)*xL)= zerosq((xL*t-1),1);
    Xzot( (t+m:sizT)*xL) = zerosq((xL*(sizT-t-m+1)),1);
    Xzot( (t:t+m-1)*xL)=X;
    T(:,t)=Xzot; 
    end
end

    case 2
     %% +++++++++++ by covariance method (T2) +++++++++++
 
 sizT=(m-n+1)*xL;
Xzot = zerosq(sizT,1);
T= zerosq(sizT,n);

for t=1:n
    Xzot=X(xL*(n-t)+1:xL*(m-t+1));
    T(:,t)   =  Xzot;     
end

    case 3
    %% +++++++++++ by covariance method (T3) +++++++++++
 
sizT=(m+1)*xL;

Xzot=zerosq(sizT,1);
T= zerosq(sizT,n);

for t=1:n
    Xzot(1:(m-n+t)*xL)=X((n-t)*xL+1:m*xL);
    Xzot((m-n+t)*xL+1:(m+1)*xL)=zerosq(xL*(n-t+1),1);
    T(:,t)   =  Xzot;     
end
 
    case 4 
 %% +++++++++++ by covariance method (T4) +++++++++++
sizT=(m+1)*xL;
Xzot=zerosq(sizT,1);
T= zerosq(sizT,n);

for t=1:n
   Xzot(t*xL+1:(m+1)*xL)=X(1:(m-t+1)*xL);
    Xzot(:,1:t*xL)=zerosq(xL*t,1); 
    T(:,t)   =  Xzot;     
end
 

end




% TT=(T'*T)/m; 
% y=[y0;zerosq(n-1,1)];
% Ty=T'*y;
