F=FILTERDATA(:,:,100:250);
NF=zeros(size(F));
%create baseline
for i=1:size(F,1)
    for j=1:size(F,2)
        PIXEL=squeeze(F(i,j,:));
        PIXEL(1:20)=PIXEL(21);%delete some pixels to create a starting point for analysis       
        %search for minimum within 30 scans after beginning of signal
        [val,idx]=min(PIXEL(1:50));
        NPIXEL=PIXEL;
        NPIXEL(1:idx)=PIXEL(21);
        NF(i,j,:)=NPIXEL;
    end
end
