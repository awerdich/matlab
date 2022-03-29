if size(ADATABASE,2)>0
    id=size(ADATABASE,2);
    size(ADATABASE)
    fprintf(['previous heart, id:',num2str(id),' name:',num2str(ADATABASE(id).name),'\n']);
    fprintf('continue? \n');
    id=input('next id:');
    pause
else
    id=1;
end
%% general information
    ADATABASE(id).age=input('age (hpf):');
    ADATABASE(id).date=input('date (yymmdd):');
    ADATABASE(id).heart=input('heart number:');
    ADATABASE(id).rate=input('HR (bpm):');
    ADATABASE(id).sinusloc=input('impulse initiation 0=a-edge 1=atrium 2=av 3=ventricle <RETURN>=fstim:');
    if isempty(ADATABASE(id).sinusloc)==1
        ADATABASE(id).stim=ADATABASE(id).rate;
    else
        ADATABASE(id).stim=input('stimulation rate 0:no stimulation (bpm):');
    end
    ADATABASE(id).comment=input('comment:','s');
    ADATABASE(id).path=datapath;
    ADATABASE(id).scanrate=scanrate;
    ADATABASE(id).name=stackfile;
    ADATABASE(id).pixelcalfactor=pixelcalfactor;%pixel calibration [um/pixel]
    ADATABASE(id).magnification=magnification;
    ADATABASE(id).DATAscanrate=scanrate/4;%scanrate of downsampled transient data
%% save tracedata
fprintf(['current id:',num2str(id),' \n']);
ADATABASE(id).name
region=input('PIXELDATA region (sa, a, av, v, vapex, vbase, o):','s');
clear regionv

if strcmp(region,'a')==1 || strcmp(region,'sa')==1 || strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
       || strcmp(region,'vbase')==1 || strcmp(region,'vapex')==1
    regionframedata=[region,'_DATA'];%ALL TRANSIENTS IN FRAME (i,j,:)
    regionmeantrace=[region,'_TRACE'];%MEAN FRAME TRANSIENT
    regionframex=[region,'_framex'];%frame width in pixels
    regionframey=[region,'_framey'];%frame height
    regionframeumwidth=[region,'_frameumwidth'];%frameumwidth
    regionframeumheight=[region,'_frameumheight'];%frameumheight
    regionframecenter=[region,'_framecenter'];%frame center locatin
    regionarea=[region,'_area'];%estimate of active tissue area in frame (SIGNALS>0)  
    regionapdlist=[region,'_apdlistms'];%list of APDs in region (ms)
    regionapd=[region,'_meanapdms',num2str(upthreshold*100)];%AVERAGE TRANSIENT DURATION (ms) 
    regionvmaxlist=[region,'_vmaxlists'];%list of VMAX values in region (1/s);
    regionvmax=[region,'_meanvmaxs'];%average maximum upstroke velocity (1/s);
    
else
    return;
end

ADATABASE(id).(regionframex)=FrameXLim;
ADATABASE(id).(regionframey)=FrameYLim;
ADATABASE(id).(regionframeumwidth)=frameumwidth;
ADATABASE(id).(regionframeumheight)=frameumheight;
ADATABASE(id).(regionframecenter)=FRAMECENTER;

%get mean transient
T=[];APDLIST=[];VMAXLIST=[];
totalscans=2000;
addscansn=totalscans;
threshold=0.5;
for i=FrameYLim(1):FrameYLim(2)
    for j=FrameXLim(1):FrameXLim(2)
        if SIGNAL(i,j)>0
            %get transient
            PIXEL=squeeze(NORMDATA(i,j,:));
            %correct baseline
            OPIXEL=PIXEL-mean(PIXEL(BASELINE));
            %re-normalize 
            [maxval,maxidx]=max(OPIXEL(UPSTROKERANGE));maxidx=maxidx+UPSTROKERANGE(1)-1;
            maxval=mean(OPIXEL(maxidx-10:maxidx+10));
            NPIXEL=OPIXEL/maxval;           
            %add data to make sure transients are long enough
            ADDSCANS=ones(addscansn,1).*mean(NPIXEL(BASELINE));
            NA=[];NA=[ADDSCANS;NPIXEL;ADDSCANS];
            maxidx=maxidx+addscansn;
            %determine upstroke and align
            k=maxidx;while NA(k)>threshold;k=k-1;end;k=k+1;%scan of >92% maximum
            beforescans=round(totalscans/4);afterscans=totalscans-beforescans-1;
            B=NA(k-beforescans:k+afterscans);
            T=[T,B];
            %get other measurements
            if APDMATRIX(i,j)>0
                APDLIST=[APDLIST;[i,j,APDMATRIX(i,j)]];
            end
            if VMAXMATRIX(i,j)>0
                VMAXLIST=[VMAXLIST;[i,j,VMAXMATRIX(i,j)]];
            end
        end
    end
end

ADATABASE(id).(regionframedata)=downsample(T,4);
ADATABASE(id).(regionmeantrace)=mean(T,2);
ADATABASE(id).(regionarea)=size(T,2);
ADATABASE(id).(regionapdlist)=APDLIST;
ADATABASE(id).(regionapd)=mean(APDLIST(:,3));
ADATABASE(id).(regionvmaxlist)=VMAXLIST;
ADATABASE(id).(regionvmax)=mean(VMAXLIST(:,3));

% save results

if size(ADATABASE,2)>1
    PREAPD=[];PREVMAX=[];TR=[];
    for i=1:size(ADATABASE,2)-1
        PREAPD=[PREAPD;ADATABASE(i).(regionapd)];
        PREVMAX=[PREVMAX;ADATABASE(i).(regionvmax)];
        TR=[TR,ADATABASE(i).(regionmeantrace)];
    end 
    figure;plot(TR,'k');hold on;plot(ADATABASE(id).(regionmeantrace),'r')
    fprintf(['mean apd (ms): ',num2str(mean(PREAPD)),' mean vmax (1/s):',num2str(mean(PREVMAX)),'\n']);
    fprintf(['new apd (ms): ',num2str(ADATABASE(id).(regionapd)),' new vmax (1/s) :',num2str(ADATABASE(id).(regionvmax)),'\n']);
end

% save to harddrive
if exist('lastname')==0
    lastname=('DATABASE');
end
newname=input(['Name of database [',num2str(lastname),']>'],'s');
if isempty(newname)==1
    filename=[lastname,'.mat'];
else
    filename=[newname,'.mat'];
end
lastname=filename(1:end-4);

input(['WRITE TO FILENAME:',num2str(filename),' [ENTER]'],'s');
save(filename,'ADATABASE');