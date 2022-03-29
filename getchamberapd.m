%determine average APD across chamber ROIs
A=[];%atrial APD
V=[];%ventricular APD
IDA=[];
IDV=[];
DATEA=[];
DATEV=[];
for id=1:length(DATABASE)
    APDA=[];
    APDA=[APDA;...
        DATABASE(id).aic_meanapd_ms20;...
        DATABASE(id).a_meanapd_ms20;...
        DATABASE(id).aoc_meanapd_ms20];
    if isempty(APDA)==0
        A=[A;max(APDA)];
        IDA=[IDA;id];
        DATEA=[DATEA;DATABASE(id).date];
    end
    
    APDV=[];
    APDV=[APDV;...
        DATABASE(id).vic_meanapd_ms20;...
        DATABASE(id).v_meanapd_ms20;...
        DATABASE(id).voc_meanapd_ms20];
    if isempty(APDV)==0
        V=[V;max(APDV)];
        IDV=[IDV;id];
        DATEV=[DATEV;DATABASE(id).date];
    end
end
        