%% The following function is preprocessing on image T.
%% T=T-mean(T)

function [Tp,meanT]=preproce(T)

[n1,n2,m]=size(T); 
sumT=zerosq(n1,n2); Tp=zerosq(n1,n2,m);

for t=1:m
    A=T(:,:,t);
    sumT=sumT+A;
end

meanT=sumT/m;

for t=1:m
    Tp(:,:,t)=T(:,:,t)-meanT;
end
    