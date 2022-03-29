if exist('DATABASE','var')==0
    DATABASE=[];id=1;
else
    id=size(DATABASE,2);
    size(DATABASE)
    fprintf(['previous heart, id:',num2str(id),' name:',num2str(DATABASE(id).name),'\n']);
    fprintf('continue? \n');
    id=input('next id:');
    pause
end
%% general information
tracedownsample=4;%keep every nth sample
if id>length(DATABASE)
    generalinput=1;
else
    DATABASE(id)
    generalinput=input('Edit general info? [RETURN]=NO, 1=YES>');
    if isempty(generalinput)==1
        generalinput=0;
    end
end

if generalinput==1
    %DATABASE(id).age=input('age (hpf):');
    DATABASE(id).name=normfile(1:end-6);
    %DATABASE(id).date=input('date (yymmdd):');
    %DATABASE(id).heart=input('heart number:');
    %DATABASE(id).rate=input('HR (bpm):');
    DATABASE(id).scanrate=scanrate;
    DATABASE(id).sinusloc=input('impulse initiation 0=atrium 1=av 2=ventricle 3=fstim:');
    %DATABASE(id).stim=input('stimulation rate 0:no stimulation (bpm):');
    %rotate heart so that the atrium is on top. OFT is either left of
    %right. This defines the surface orientation of the heart
    %DATABASE(id).orientation=input('With A on top, OFT is LEFT(1) or RIGHT(2):');
    %DATABASE(id).comment=input('comment:','s');
    DATABASE(id).path=datapath;
    DATABASE(id).pixelcalfactori=pixelcalfactor_i;%pixel calibration [um/pixel]
    DATABASE(id).pixelcalfactorj=pixelcalfactor_j;%pixel calibration [um/pixel]
    DATABASE(id).magnification=magnification;
    DATABASE(id).sideport=sideport;
    DATABASE(id).binning=binning;
    DATABASE(id).DATAscanrate=scanrate/tracedownsample;%scanrate of downsampled transient data
    DATABASE(id).filterlp=lpcutoff;%low pass cutoff frequency
    DATABASE(id).median=medianfilter;%median filter degree
    DATABASE(id).softwareversion=pwd;%software version
    if exist('golaywindow','var')==1
        DATABASE(id).golaywindowAPDMAP=golaywindow;
    end
    DATABASE(id).filterspatpix=spatpix;%spatial filter matrix
    if exist('angleunitvector','var')==1
        DATABASE(id).ANGLEUNITVECTORXY=angleunitvector;%angle unit vector for specific region
    end
    if exist('POL','var')==1
        DATABASE(id).pol1repol2=pol1repol2;%1=POL 2=REPOL
        %save POL and REPOL matrices
        DATABASE(id).POL=POL;
        DATABASE(id).REPOL=REPOL;
    end
    DATABASE(id).SIGNAL=SIGNAL;
    %add stimulus field
    if exist('EPISODESTIM','var')==1 && isempty(EPISODESTIM)==0
        DATABASE(id).STIM=[];
    end
end
%% save propagation velocity data

fprintf(['current id:',num2str(id),' \n']);
DATABASE(id).name
fprintf('continue? \n');
region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
clear regionv


if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
|| strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
|| strcmp(region,'vic')==1 || strcmp(region,'voc')==1
   %FIELD DEFINITIONS
   regionmaxrmse=[region,'_maxrmse'];%maxRMSE
   regionvelinterval=[region,'_VELINTERVAL'];%velocity inteval for fitted vectors
   regionv=[region,'_vel'];%mean propagation velocity
   regionvstd=[region,'_velstd'];%standard deviation of FINAL vector list (SORTV)
   regionsignalframe=[region,'_SIGNALFRAME'];
   regionframearea=[region,'_framearea_um'];%active tissue area in frame (um^2, from SIGNALFRAME)
   regiontotalarea=[region,'_totalarea_um'];%estimated total tissue area (um^2, from SIGNAL)
   regionfits=[region,'_numfits'];%total number fits performed in frame
   regionwindowx=[region,'_windowx'];%frameumwidth
   regionwindowy=[region,'_windowy'];%frameumheight
   regionwindowcenter=[region,'_FRAMECENTER'];%frame center location
   regionframevi=[region,'_VFRAMEi'];%conduction velocity vector i-component [pixel/frame]
   regionframevj=[region,'_VFRAMEj'];%conduction velocity vector j-component [pixel/frame]
   regionframex=[region,'_XFRAMEINTERP0'];%x-coordinates of frame region
   regionframey=[region,'_YFRAMEINTERP0'];%y-coordinates of frame region
   regionsignalframeinterp=[region,'_SIGNALFRAMEINTERP0'];%Binary image: Signals in frame
   regionplotbin=[region,'_PLOTVELOCITYBIN'];%binary image, 1:pixel used 0:pixel filtered
   regionsignalframeinterpolygon=[region,'_POLYROI'];%Binary image: Signals in polygon
   regionrmse=[region,'_RMSE'];%RMSE for each vector in frame
   regionpol=[region,'_POL'];
   regionrepol=[region,'_REPOL'];
end

DATABASE(id).(regionplotbin)=VPLOTBIN;%binary image to mark pixels used for velocity calculation
DATABASE(id).(regionmaxrmse)=maxrmse;%marmse
DATABASE(id).(regionvelinterval)=VELINTERVAL;
DATABASE(id).(regionv)=meanvelocity;%mean propagation velocity
DATABASE(id).(regionvstd)=stdvelocity;%standard deviation of FINAL vector list (SORTV)
DATABASE(id).(regionsignalframe)=SIGNALFRAME;
DATABASE(id).(regionframearea)=framearea;%active tissue area in frame (um^2, from SIGNALFRAME) 
DATABASE(id).(regiontotalarea)=area;%estimated total tissue area (um^2, from SIGNAL)
DATABASE(id).(regionfits)=framepixels;%total number fits performed in frame
DATABASE(id).(regionwindowx)=frameumwidth;
DATABASE(id).(regionwindowy)=frameumheight;
DATABASE(id).(regionwindowcenter)=FRAMECENTER;
DATABASE(id).(regionframevi)=VFRAMEi;
DATABASE(id).(regionframevj)=VFRAMEj;
DATABASE(id).(regionframex)=XFRAMEINTERP_0;
DATABASE(id).(regionframey)=YFRAMEINTERP_0;
DATABASE(id).(regionsignalframeinterp)=SIGNALFRAMEINTERP_0;
if exist('POLYROI','var')==1
    DATABASE(id).(regionsignalframeinterpolygon)=POLYROI;
end
DATABASE(id).(regionrmse)=RMSEFRAME;

%extract POL and REPOL times in ROIs
%get the coordinates and initialize data ROIs
ROI=SIGNAL;
X=[FrameXLim(1):FrameXLim(2)]';
Y=[FrameYLim(1):FrameYLim(2)]';
POL_ROI=zeros(length(Y),length(X));
REPOL_ROI=zeros(size(POL_ROI));
for i=Y(1):Y(end)
    k=i-Y(1)+1;
    for j=X(1):X(end)
        l=j-X(1)+1;
        if SIGNAL(i,j)==1
            ROI(i,j)=0;
            POL_ROI(k,l)=POL(i,j);
            REPOL_ROI(k,l)=REPOL(i,j);
        end
    end
end
%save data in DATABASE
DATABASE(id).(regionpol)=POL_ROI;
DATABASE(id).(regionrepol)=REPOL_ROI;


%update angleunitvector if it has been changed in showcontours
if isfield(DATABASE,'ANGLEUNITVECTORXY')==1 && (isempty(DATABASE(id).ANGLEUNITVECTORXY)==1 || sum(DATABASE(id).ANGLEUNITVECTORXY==angleunitvector)<2)
   DATABASE(id).ANGLEUNITVECTORXY=angleunitvector;%angle unit vector for specific region
   fprintf('Unitvector changed.\n');pause
end


%display data
%DATABASE(id);
%collect previous data
V=[];ASTD=[];
for i=1:length(DATABASE)
    if DATABASE(i).sinusloc==DATABASE(id).sinusloc
        V=[V;DATABASE(i).(regionv)];
    end
end
V
fprintf(['mean velocity for previous data: ',num2str(mean(V)),' new value: ',num2str(DATABASE(id).(regionv)),' \n']);
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

vinput=input(['WRITE TO FILENAME:',num2str(filename),' [ENTER]'],'s');
%save(filename,'DATABASE');


%% SAVE APDs
if exist('APDMATRIX')==1
    if exist('vinput','var')==0
        region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
    end
    if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
        || strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
        || strcmp(region,'vic')==1 || strcmp(region,'voc')==1
        %NEW FIELD DEFINITIONS
        regionframedata=[region,'_TRACEDATA'];%ALL T TRANSIENT
        regionmeantrace=[region,'_TRACE'];%MEAN FRAME TRANSIENT
        regionapdlist=[region,'_APDLIST_ms'];%list of APDs in region (ms)
        regionapd=[region,'_meanapd_ms',num2str(apdthreshold*100)];%AVERAGE TRANSIENT DURATION (ms) 
        regionvmaxlist=[region,'_VMAXLIST_s'];%list of VMAX values in region (1/s);
        regionvmaxlistsm=[region,'_VMAXLISTSM_s'];%list of VMAX values from smoothed transients (1/s);
        regionvmax=[region,'_meanvmax_s'];%average maximum upstroke velocity (1/s);
        regionvmaxsm=[region,'_meanvmaxsm_s'];%average vmax from smoothed transients(1/s)
        regionstimtrace=[region,'_STIMTRACE'];%stimulus transient;
        
        %some information that needs to be saved if CVs are not calculated
        regionwindowx=[region,'_windowx'];%frameumwidth
        regionwindowy=[region,'_windowy'];%frameumheight
        regionwindowcenter=[region,'_FRAMECENTER'];%frame center location
        
        DATABASE(id).(regionwindowx)=frameumwidth;
        DATABASE(id).(regionwindowy)=frameumheight;
        DATABASE(id).(regionwindowcenter)=FRAMECENTER;
    else
        fprintf('Wrong input!\n')
        return
    end

%get mean transient
T=[];APDLIST=[];VMAXLIST=[];VMAXLISTSM=[];KLIST=[];
totalscans=5000;
addscansn=totalscans;
threshold=0.5;
for i=FrameYLim(1):FrameYLim(2)
    for j=FrameXLim(1):FrameXLim(2)
        if SIGNAL(i,j)>0
            %get transient
            PIXEL=squeeze(NORMDATA(i,j,:));
            BASELINE=squeeze(BASELINESCANS(i,j,:));
            %correct baseline
            OPIXEL=PIXEL-mean(PIXEL(BASELINE));
            %re-normalize 
            [maxval,maxidx]=max(OPIXEL(UPSTROKERANGE));maxidx=maxidx+UPSTROKERANGE(1)-1;
            maxval=mean(OPIXEL(maxidx-10:maxidx+10));
            NPIXEL=OPIXEL/maxval;           
            %add data to make sure transients are long enough
            ADDSCANS=ones(addscansn,1).*mean(NPIXEL(BASELINE));
            NA=[];NA=[ADDSCANS;NPIXEL;ADDSCANS];
            if exist('EPISODESTIM','var')==1 && isempty(EPISODESTIM)==0 && i==FrameYLim(1) && j==FrameXLim(1)
                NSTIM=[ADDSCANS;EPISODESTIM;ADDSCANS];
            end
            maxidx=maxidx+addscansn;
            %determine upstroke and align
            k=maxidx;while NA(k)>threshold;k=k-1;end;k=k+1;%scan of >92% maximum
            beforescans=floor(round(totalscans/4));afterscans=totalscans-beforescans-1;
            B=NA(k-beforescans:k+afterscans);
            T=[T,B];
            KLIST=[KLIST;k];%save index of maximum to determine mean index in ROI
            %get other measurements from saved matrices
            if APDMATRIX(i,j)>0
                APDLIST=[APDLIST;[i,j,APDMATRIX(i,j)]];
                VMAXLIST=[VMAXLIST;[i,j,VMAXMATRIX(i,j)]];
            end
            if exist('VMAXMATRIXSM','var')==1 && VMAXMATRIXSM(i,j)>0
                VMAXLISTSM=[VMAXLISTSM;[i,j,VMAXMATRIXSM(i,j)]];
            end
       end
    end
end

%add stimulus if available
if exist('NSTIM','var')==1
    if isempty(DATABASE(id).STIM)==1
        DATABASE(id).STIM=NSTIM(round(mean(KLIST))-beforescans:round(mean(KLIST))+afterscans);
    end
end

%DATABASE(id).DATAscanrate=scanrate/tracedownsample;
tracedownsample=DATABASE(id).scanrate/DATABASE(id).DATAscanrate;

DATABASE(id).(regionframedata)=downsample(T,tracedownsample);
DATABASE(id).(regionmeantrace)=mean(T,2);
DATABASE(id).(regionapdlist)=APDLIST;
DATABASE(id).(regionapd)=mean(APDLIST(:,3));
DATABASE(id).(regionvmaxlist)=VMAXLIST;
DATABASE(id).(regionvmaxlistsm)=VMAXLISTSM;
DATABASE(id).(regionvmax)=mean(VMAXLIST(:,3));
if length(VMAXLISTSM)>1
    DATABASE(id).(regionvmaxsm)=mean(VMAXLISTSM(:,3));
else
    DATABASE(id).(regionvmaxsm)=0;
end

%save stimulus


% lookup previous measurements
PREAPD=[];PREVMAX=[];PREVMAXSM=[];TR=[];
if size(DATABASE,2)>1
    for i=1:size(DATABASE,2)
        if DATABASE(i).sinusloc==DATABASE(id).sinusloc
            PREAPD=[PREAPD;DATABASE(i).(regionapd)];
            PREVMAX=[PREVMAX;DATABASE(i).(regionvmax)];
            PREVMAXSM=[PREVMAXSM;DATABASE(i).(regionvmaxsm)];
            TR=[TR,DATABASE(i).(regionmeantrace)];
        end
    end
    figure;plot(TR,'k');hold on;
    fprintf(['Mean previous APD (ms): ',num2str(mean(PREAPD)),'\n']);
    fprintf(['Mean Vmax from smoothed transients (1/s): ',num2str(mean(PREVMAXSM)),'\n']);
end
if isempty(TR)==1
    figure
end
plot(DATABASE(id).(regionmeantrace),'r')

%plot upstroke velocity for FRAMECENTER
i=round(FRAMECENTER(1));j=round(FRAMECENTER(2));
FCENTERTRANS=squeeze(NORMDATA(i,j,:));
if exist('BASELINESCANS','var')==1
    BASELINE=squeeze(BASELINESCANS(i,j,:));
end
[upidx,downidx,T]=apd(FCENTERTRANS,apdthreshold,UPSTROKERANGE,BASELINE,scanrate);
if exist('golaywindow','var')==0
    golaywindow=31;%smoothing filter window sent to CAM on 10/6/2010 for cplexapd
end
[D0,D1,D2]=smoothdiff(T(UPSTROKERANGE),golaywindow);
figure;plot(D1*scanrate);
PREAPD
fprintf([' New APD (ms): ',num2str(DATABASE(id).(regionapd)),'\n']);
fprintf([' New Vmax (1/s): ',num2str(DATABASE(id).(regionvmaxsm)),'\n']);



% save to harddrive
if exist('lastname')==0
    lastname=('DATABASE');
end
%newname=input(['Name of database [',num2str(lastname),']>'],'s');
if isempty(newname)==1
    filename=[lastname,'.mat'];
else
    filename=[newname,'.mat'];
end
lastname=filename(1:end-4);

%input(['WRITE TO FILENAME:',num2str(filename),' [ENTER]'],'s');
save(filename,'DATABASE'); 
end
DATABASE(id)