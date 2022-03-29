function NDATA=ndata(DATA,BASELINESCANS)
NDATA=zeros(size(DATA));
%normalize DATA using BASELINESCANS
for i=1:size(DATA,2)
    T=DATA(:,i);
    OT=T-mean(T(BASELINESCANS));
    NDATA(:,i)=OT/max(OT(180:250));
end
