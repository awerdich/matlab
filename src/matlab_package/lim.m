%determine transient limits
DATA=[CONTROL(:,1:4),CYANIDE(:,1:4)];
naverage=20;
%get all amplitudes
for i=1:size(DATA,2)
    ALLAMP=DATA(:,i);
    ALLAMP=[[1:length(ALLAMP)]',ALLAMP];
    %sort amplitudes
    SORTAMP=sortrows(ALLAMP,2);
    %determine mean extrema
    MAXAMP(i)=mean(SORTAMP(end-naverage+1:end,2));
    MINAMP(i)=mean(SORTAMP(1:naverage,2));
end