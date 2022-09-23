ROITIME=[];
avscans=5;
for id=1:length(DATABASE)
    if isempty(DATABASE(id).voc_POL)==0
        scanrate=DATABASE(id).scanrate;
        %activation matrix
        POL=DATABASE(id).voc_POL;
        %collect activation times
        T=[];
        for i=1:size(POL,1)
            for j=1:size(POL,2)
                if POL(i,j)>0
                    T=[T;POL(i,j)];
                end
            end
        end
        %sort T
        ST=sortrows(T,1);
        %determine activation time
        ROITIME=[ROITIME;[id,DATABASE(id).rate,(mean(ST(end-avscans+1:end,1))-mean(ST(1:avscans,1)))/scanrate*1.0e3]];
    end
end
save ROITIME.txt ROITIME -ascii -tabs