
function DM=DistanceMatrix(dsites,ctrs)
[M,s]=size(dsites);[N,s]=size(ctrs);
DM=zeros(M,N);
%Accumulate sum of squares of coordinate disserences
%The ndgrid command produces two MxN matrices;

for d=1:s
    [dr cc]=ndgrid(dsites(:,d),ctrs(:,d));
    DM=DM+(dr-cc).^2;
end
DM=sqrt(DM);