function T=QuaCirculant(col)
m=size(col,1);
 CM = zerosq(m,m); %Initialize Circulant Matrix (CM)
CM(:,1) =col; %place the input values in the first column of CM
for kk = 1:(size(CM,1) - 1)
    CM(1,kk+1) = CM(end,kk);
    CM(2:end,kk+1) = CM( 1:end - 1,kk);
end
T=CM;
end

