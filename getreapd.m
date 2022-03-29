%get APDs from DATABASE at new level
%% load database
[file,path]=uigetfile('*.mat','Load APD database');
load([path,file]);
% define parameters
threshold=0.5;
BASELINE=[150:250];
UPSTROKERANGE=[250:450];
scanrate=DATABASE(1).DATAscanrate;
REGS={'a','aic','aoc','av','v','vic','voc','o'};
A=[];AIC=[];AOC=[];AV=[];V=[];VIC=[];VOC=[];OFT=[];
 hdl = waitbar(0,'PROCESSING DATABASE ENTRIES');
for id=1:length(DATABASE)
for reg=1:length(REGS)
    %transient array for region
    fieldname=[REGS{reg},'_TRACEDATA'];
    if isfield(DATABASE,fieldname)==1 && isempty(DATABASE(id).(fieldname))==0
        TRACEDATA=DATABASE(id).(fieldname);
        if isempty(TRACEDATA)==0
            %run apd function on this
            [UPIDX,DOWNIDX,MEANTRANS]=apd(TRACEDATA,threshold,UPSTROKERANGE,BASELINE,scanrate);
            %calculate APDs
            APD=mean((DOWNIDX-UPIDX)/scanrate*1.0e3);
            %put the data into the correct VECTOR
            if reg==1
                A=[A;APD;];
            elseif reg==2
                AIC=[AIC;APD];
            elseif reg==3
                AOC=[AOC;APD];
            elseif reg==4
                AV=[AV;APD];
            elseif reg==5
                V=[V;APD];
            elseif reg==6
                VIC=[VIC;APD];
            elseif reg==7
                VOC=[VOC;APD];
            elseif reg==8
                OFT=[OFT;APD];
            end
        end
    end
end
waitbar(id/length(DATABASE))
end
close(hdl)
%% save results
folder=uigetdir(path,'Save APD text files');
save([folder,'\A.txt'],'A','-ascii')
save([folder,'\AIC.txt'],'AIC','-ascii');
save([folder,'\AOC.txt'],'AOC','-ascii');
save([folder,'\AV.txt'],'AV','-ascii');
save([folder,'\V.txt'],'V','-ascii');
save([folder,'\VIC.txt'],'VIC','-ascii');
save([folder,'\VOC.txt'],'VOC','-ascii');
save([folder,'\OFT.txt'],'OFT','-ascii');