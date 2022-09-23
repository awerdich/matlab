%% calculate activation times by searching pixel timeseries
% requires NORMDATA, SIGNALPIXELS
DATA=NORMDATA;
SIGNALS=ASIGNALPIXELS;
upthresh=0.5;%depolarization threshold
downthresh=0.2;%repolarization threshold
UFITINTERVAL=[upthresh-0.2;upthresh+0.2];%fit interval in upstroke
DFITINTERVAL=[downthresh-0.1;downthresh+0.3];%fit interval in decay
POL=zeros(size(DATA,1),size(DATA,2));%time of activation
REPOL=zeros(size(POL));%time of repolarization
ACTIVATIONTIMEFAIL=zeros(size(POL));%marks pixels that failed to produce an activation time

if exist ('BASELINE')==0 || exist('UPSTROKERANGE')==0
    if exist('FBR')==0
        FBR=[];
    end
    [PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,SCATTER,ASIGNALPIXELS,DATA);
    % select range based on sample pixels (SPIXELS)
    figure('name','range of upstrokes (2 clicks), baseline (2 clicks)')
    SPIXELS=[squeeze(DATA(PIXELCOORDINATES(1,2),PIXELCOORDINATES(1,1),:)),squeeze(DATA(PIXELCOORDINATES(2,2),PIXELCOORDINATES(2,1),:))];
    plot(SPIXELS(:,1),'k');hold on
    plot(SPIXELS(:,2),'r');
    fprintf('select range of upstrokes (2 clicks)\n');
    [ax,ay]=ginput(2);
    UPSTROKERANGE=[round(ax(1)):round(ax(2))]';
    plot(UPSTROKERANGE,SPIXELS(UPSTROKERANGE,:),'m');
    fprintf('select baseline (2 clicks)\n');
    [ax,ay]=ginput(2);
    BASELINE=[round(ax(1)):round(ax(2))]';
    plot(BASELINE,SPIXELS(BASELINE,:),'g');
end
%% search loop
tstart=tic;%start timer
hdl = waitbar(0,['searching for upstrokes and repolphases']);
for i=1:size(DATA,1)
    for j=1:size(DATA,2)    
%% get pixel time series
        if SIGNALS(i,j)==1 
            PIXEL=squeeze(DATA(i,j,:));
            %correct baseline
            BASELINE=squeeze(BASELINESCANS(i,j,:));
            OPIXEL=PIXEL-mean(PIXEL(BASELINE));
            %re-normalize 
            [maxval,maxidx]=max(OPIXEL(UPSTROKERANGE));maxidx=maxidx+UPSTROKERANGE(1)-1;
            NPIXEL=OPIXEL/mean(OPIXEL(maxidx-1:maxidx+3));
            
            %start searching from maximum
            up=maxidx;while NPIXEL(up)>upthresh && up>1;up=up-1;end
            %find fit interval boundaries: start with low border of fit interval
            upfitlow=up;while NPIXEL(upfitlow)>UFITINTERVAL(1) && upfitlow>1;upfitlow=upfitlow-1;end
            if upfitlow>1
                upfithigh=up;while NPIXEL(upfithigh)<UFITINTERVAL(2);upfithigh=upfithigh+1;end
                UPXDATA=[upfitlow:upfithigh]';UPYDATA=NPIXEL(UPXDATA);
                %fit
                if length(UPXDATA)>5
                    fitorder=5;
                else
                    fitorder=1;
                end
                [upfit,S,MU]=polyfit(UPXDATA,UPYDATA,fitorder);
                %calculate threshold index from fitted polynom in finer
                UPXDATAFIT=[UPXDATA(1):0.001:UPXDATA(end)]';
                UF=(UPXDATAFIT-MU(1))./MU(2);%coordinate transformation to calculate fitted polynom
                UPYDATAFIT=polyval(upfit,UF);
                
                %search for new threshold on fit
                upfitidx=1;while UPYDATAFIT(upfitidx)<upthresh && upfitidx<length(UPYDATAFIT);upfitidx=upfitidx+1;end
                upidx=UPXDATAFIT(upfitidx);
                POL(i,j)=upidx;
            end
            
            
            if POL(i,j)>0 %if you cannot fit the upstroke, the decay wont work either!
                %filter transient [FDATA]=fdtrans(RDATA,lpcutoff,medscans,scanrate)
                FPIXEL=fdtrans(NPIXEL,200,5,scanrate);
                down=maxidx;while FPIXEL(down)>downthresh && down<length(FPIXEL);down=down+1;end
                %find fit interval boundaries: start with high border
                downfithigh=down;while FPIXEL(downfithigh)<DFITINTERVAL(2) && downfithigh>maxidx;downfithigh=downfithigh-1;end
                downfitlow=down;while FPIXEL(downfitlow)>DFITINTERVAL(1) && downfitlow<length(FPIXEL);downfitlow=downfitlow+1;end
                DOWNXDATA=[downfithigh:downfitlow]';DOWNYDATA=FPIXEL(DOWNXDATA);
                if length(DOWNXDATA)>3
                    [downfit,S,MU]=polyfit(DOWNXDATA,DOWNYDATA,3);
                    %threshold index from fitted polynom
                    DOWNXDATAFIT=[DOWNXDATA(1):0.001:DOWNXDATA(end)]';            
                    DF=(DOWNXDATAFIT-MU(1))./MU(2);
                    DOWNYDATAFIT=polyval(downfit,DF);
                    %search for new index on fit
                    downfitidx=1;while DOWNYDATAFIT(downfitidx)>downthresh & downfitidx<length(DOWNYDATAFIT);downfitidx=downfitidx+1;end
                    downidx=DOWNXDATAFIT(downfitidx);
                else
                    downidx=0;
                end
            else
                downidx=0;
            end
            REPOL(i,j)=downidx;
          end
        end
%% Start over with next pixel
    waitbar(i/size(DATA,1));
end
close(hdl)
elapsedtime=toc(tstart);fprintf(['Elapsed time: ',num2str(elapsedtime),' s\n']);
%% generate ACTIVATIONLIST
LIST=[];ACTIVATIONLIST=[];
for i=1:size(POL,1)
    for j=1:size(POL,2)
        if SIGNALS(i,j)==1 && SCATTER(i,j)==1
            LIST=[LIST;[i,j,POL(i,j)]];
        end
    end
end
%sort activationlist
ACTIVATIONLIST=zeros(size(LIST));
[~,IX]=sort(LIST(:,3),'ascend');
for i=1:length(IX)
    ACTIVATIONLIST(i,:)=LIST(IX(i),:);
end
%% determine first and last pixels that are inside plotting area and unitvector
%LIST=[LIST;[i,j,POL(i,j)]];
SIGNALS(SCATTER==0)=0;
%remove border pixels from SIGNALS
B=bwperim(SIGNALS);
BSIGNALS=SIGNALS;BSIGNALS(B==1)=0;

k=1;i=ACTIVATIONLIST(k,1);j=ACTIVATIONLIST(k,2);
while BSIGNALS(i,j)==0 && k<50 
    k=k+1;
    i=ACTIVATIONLIST(k,1);j=ACTIVATIONLIST(k,2);
end
first=k; 
k=size(ACTIVATIONLIST,1);i=ACTIVATIONLIST(k,1);j=ACTIVATIONLIST(k,2);
while BSIGNALS(i,j)==0 && k>size(ACTIVATIONLIST,1)-50
    k=k-1;
    i=ACTIVATIONLIST(k,1);j=ACTIVATIONLIST(k,2);
end
last=k;
%% set activation times for pixels that could not be measured
SIGNAL=ASIGNALPIXELS;
SIGNAL(SCATTER==0)=0;
ACTIVATIONTIMEFAIL(SIGNAL==0)=0;
FIXPOL=POL;%final fixed activation matrix
C=POL;%temporary matrix
WIENC=wiener2(C,[6,6]);
MEDC=medfilt2(WIENC,[6,6]);
for i=1:size(ACTIVATIONTIMEFAIL,1)
    for j=1:size(ACTIVATIONTIMEFAIL,2)
        if ACTIVATIONTIMEFAIL(i,j)>0
            FIXPOL(i,j)=MEDC(i,j);
        end
    end
end
%% save avtivation times
polfile=[normfile(1:end-5),'ATIME.mat'];
if exist('datapath','var')==0 || exist(datapath,'file')==0
    [file datapath]=uiputfile('*.*','SELECT DATA FOLDER');
end
save([datapath,polfile],'POL','REPOL','ACTIVATIONLIST','upthresh',...
    'downthresh','BASELINE','UPSTROKERANGE','datapath','first','last');