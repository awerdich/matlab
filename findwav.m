% return wavelengths in file and scan intervals
function [WAVELENGTHS,SCANINTERVALS]=findwav(WAV,READY)
%% find wavelength intervals

if isempty(READY)==1
    %find plateau in WAV if READY signal not available to determine
    %wavelength intervals based on 3-point plateaus
    wavtolerance=30;%2nm wavelength uncertainty within constant wavelength
    WAVIDX=[];
    i=2;
    while i<length(WAV)
        if abs(WAV(i-1)-WAV(i))<wavtolerance && abs(WAV(i+1)-WAV(i))<wavtolerance
            %plateau with three points
            j=i+1;
            while abs(WAV(j)-WAV(i))<wavtolerance && j<length(WAV)-1,j=j+1;end
            WAVIDX=[WAVIDX;[i+1,j-1,mean(WAV(i+1:j-1))]];
            i=j;%next search starts after plateau
        else
            i=i+1;
        end
    end

elseif isempty(READY)==0;%READY SIGNAL AVAILABLE; use it to find wavelength intervals
    WAVIDX=[];
    readythreshold=4250;
    i=2;
    while i<length(READY)-1
        if READY(i)>readythreshold
            %found READY plateau, determine length
            left=i;while READY(i)>readythreshold && i<length(READY)-1,i=i+1;end;right=i-1;
            if READY(right)>readythreshold && (right-left)>=1 && i<length(READY)-5%last index could be due to end of file (EOF)
                WAVIDX=[WAVIDX;[right-1,right-1,WAV(right)]];
                %WAVIDX=[WAVIDX;[left+1,right,mean(WAV(left:right))]];
            end
        end
        i=i+1;
    end
end
   
%% find wavelengths
LAMBDA=round(DAC(WAVIDX(1,3)));%WAVELENGTHS [nm]
dl=7;%[nm] two wavelengths are considered different if their values differ by >dl
for i=2:size(WAVIDX,1)
    l=round(DAC(WAVIDX(i,3)));%wavelength under consideration
    if l-LAMBDA(end)>dl%wavelength is larger than last in list
        LAMBDA=[LAMBDA;l];%add new wavelength to end of list
    elseif LAMBDA(1)-l>dl%wavelength is smaller than first in list
        LAMBDA=[l;LAMBDA];%insert new wavelength to beginning of list
    end
end
%% sort wavelength time intervals
LAMBDAIDX=zeros(length(WAVIDX),2*length(LAMBDA));
for i=1:size(WAVIDX,1)-1
    assigned=0;
    for j=1:length(LAMBDA)
        if abs(DAC(WAVIDX(i,3))-LAMBDA(j))<dl
            lambdaidxfirstcolumn=2*j-1;%first column in LAMBDAIDX 
            assigned=1;%wavelength in WAVIDX was assigned to wavlength in LAMBDA
        end
        if j==length(LAMBDA) && assigned==0
            sprintf(['wavelength on WAVIDX(',num2str(i),',3) could not be assigned\n'])
            return
        end 
    end


    %determine last row in LAMBDAIDX and enter corresponding time intervals
    %given the interval in WAVIDX and the lambdafirstcolumn, enter scan
    %interval into LAMBDAIDX database
    INTERVAL=[WAVIDX(i,1);WAVIDX(i,2)];
        
    if LAMBDAIDX(1,lambdaidxfirstcolumn)==0
        LAMBDAIDX(1,lambdaidxfirstcolumn:lambdaidxfirstcolumn+1)=[INTERVAL(1),INTERVAL(2)];
    else
        row=1;while LAMBDAIDX(row,lambdaidxfirstcolumn)~=0,row=row+1;end
        LAMBDAIDX(row,lambdaidxfirstcolumn:lambdaidxfirstcolumn+1)=[INTERVAL(1),INTERVAL(2)];
    end
    LASTROW=[];
    for k=1:2:size(LAMBDAIDX,2)
        m=1;while LAMBDAIDX(m,k)>0;m=m+1;end;LASTROW=[LASTROW;m-1];
    end
end
%% plot results selection
figure
PI=[600:800];
%plot wavelength signal
plot(PI,WAV(PI));set(gca,'XLim',[PI(1),PI(end)]);hold on;plot(PI,WAV(PI),'r+');
if isempty(READY)==0
    %offset and normalize READY DATA
    RD=(((READY-min(READY))/max(READY))*(max(WAV)-min(WAV)))+min(WAV);
    %plot READY SIGNAL
    plot(PI,RD(PI),'k');set(gca,'XLim',[PI(1),PI(end)]);hold on;plot(PI,RD(PI),'k+');
end

    
for i=1:size(WAVIDX,1)
    if PI(1)<=WAVIDX(i,1) && WAVIDX(i,2)<=PI(end)
        plot(WAVIDX(i,1):0.1:WAVIDX(i,2),WAVIDX(i,3),'k.','MarkerSize',5)
    end
end
i=1;
%find first index in LAMBDAIDX larger than PI(1)
while min(LAMBDAIDX(i,:))<PI(1),i=i+1;end;i=i-1;
while max(LAMBDAIDX(i,:))<=PI(end)
    %plot intervals for first wavelength
    k=LAMBDAIDX(i,1);l=LAMBDAIDX(i,2);
    plot(k:l,WAV(k:l),'ko')
    %plot intervals for second wavelength
    m=LAMBDAIDX(i,3);n=LAMBDAIDX(i,4);
    plot(m:n,WAV(m:n),'mo')
    i=i+1;
end
%% return results
%cut LAMBDAIDX
i=1;
while min(LAMBDAIDX(i,:))>0,i=i+1;end,lastrow=i-2;   
WAVELENGTHS=LAMBDA;
SCANINTERVALS=LAMBDAIDX(1:lastrow,:);