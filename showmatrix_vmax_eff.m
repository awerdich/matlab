%autoset mapconfig to 20,1,0.5,enter, for efficiency
mapconfig_eff_voltage

description=input('DESCRIPTION>','s');

%% user input
if exist('pixelcalfactor_x')==0
    mapconfig
end
    
clear FRAMECENTER;%close all
%DIASTOLIC CA
%C=RLOW;limits=1;

%CA RELEASE
%C=(RHIGH-RLOW);limits=3;

%CA TRANSIENT DURATION                             
%C=DUR;limits=4;

%ACTION POTENTIAL DURATIONS
%C=APDMATRIX;limits=5;

%ACTION POTENTIAL UPSTROKE VELOCITIES
C=VMAXMATRIXSM;limits=6;

%CONDUCTION VELOCITIES
%C=VELMATRIX;limits=7;

%ANGLES
%C=ANG180;limits=8;

%S2 INTERVALS
%C=S2AMPMATRIXBIN;limits=9;

%DIVERGENCE
%C=cvdivergence(VFRAMEj,VFRAMEi);limits=10;
%ASIGNALPIXELS = SIGNALFRAMEINTERP_0;

showframe=1;%1=show frame in isochronal map

setframeumwidth=30;%frame unit length [um]. If zero, click corners of frame
setframeumheight=30;

opengl software;
SPATWIEN=[5 5];
SPATMED=[5 5];
delborder=3;%delete pixels at border to eliminate filtering artifacts
outputresolution=[1024,1024];%output image resolution [height length]

if exist('magnification','var')==0;mapconfig;end
pixelcalfactor=pixelcalfactor_x;
unitlength=20;%length of scalebar in um

showcolorbar=0;
showscalebar=0;

%COLORBAR POSITION
lright=[0.75 0.02 0.05 0.3];%right side colorbar 'position'
lleft=[0.01 0.02 0.05 0.3];%left side colorbar 'position'
uleft=[0.01 0.68 0.05 0.3];
uright=[0.72 0.68 0.05 0.3];
colorbarposition=uleft;

%R LIMITS

%LIMITS LOW
if limits==1
    %STANDARD
    lowlimit=0.5;
    highlimit=2.5;
    stepsize=0.1;
    file=[stackfile(1:end-4),'-LOW',num2str(lowlimit),'-',num2str(highlimit)];
    %lowlimit=0.8;
    %highlimit=1.2;
    %stepsize=0.1;

%LIMITS HIGH
elseif limits==2
    lowlimit=1.5;
    highlimit=2.5;
    stepsize=0.1;
    file=[stackfile(1:end-4),'-HIGH',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS RELEASE
elseif limits==3
    lowlimit=0.0;
    highlimit=0.10;
    stepsize=0.01;
    file=[stackfile(1:end-4),'-REL',num2str(lowlimit),'-',num2str(highlimit)];
    
%LIMITS CA DURATION
elseif limits==4
    lowlimit=150;
    highlimit=500;
    stepsize=50;
    file=[stackfile(1:end-4),'-DUR',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS APD
elseif limits==5;
    lowlimit=50;
    highlimit=600;
    stepsize=50;
    file=[stackfile(1:end-4),'-APD',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS VMAX
elseif limits==6
    lowlimit=20;
    highlimit=120;
    stepsize=20;
    file=[stackfile(1:end-4),'-VMAX',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS CONDUCTION VELOCIIES
elseif limits==7
    lowlimit=0;
    highlimit=50;
    stepsize=5;
    file=[stackfile(1:end-4),'-CV',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS ANGLES
elseif limits==8
    lowlimit=0;
    highlimit=180;
    stepsize=90;
    file=[stackfile(1:end-4),'-ANG',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS S2AMPLITUDES
elseif limits==9
    lowlimit=0.1;
    highlimit=1.0;
    stepsize=0.2;
    file=[stackfile(1:end-4),'-S2AMP',num2str(lowlimit),'-',num2str(highlimit)];

%LIMITS DIVERGENCE
elseif limits==10
    lowlimit=-1;
    highlimit=+1;
    stepsize=0.2;
    file=[stackfile(1:end-4),'-DIV',num2str(lowlimit),'-',num2str(highlimit)];
end

%resizefactor
if exist('SHIFTB')==1
    resizefactor=fresizefactor;
else
    resizefactor=10;
end

%frame
%NFRAME=[];%comment to save frame coordinates
sframe=0;
framecolor=[1 0 0];%color of frame
frameweight=2;
showcenter=0;%1=show only frame center
%% filter data matrix
if sum(SPATWIEN)==0 || isempty(SPATWIEN)==1
    WIENC=C;
else
    WIENC=wiener2(C,SPATWIEN);
end

if sum(SPATMED)==0 || isempty(SPATMED)==1
    MEDC=WIENC;
else
    MEDC=medfilt2(WIENC,SPATMED);
end
%% replace border of active pixels 
%combine ASIGNALPIXELS and SCATTER in SIGNAL
 %SIGNAL MATRIX FOR WHOLE FIELD
SIGNAL=ASIGNALPIXELS;
SIGNAL(SCATTER==0)=0;
SIGNALCOPY=SIGNAL;

BORDER=zeros(size(SIGNAL));
if delborder>0
    for j=1:delborder
        B=bwboundaries(SIGNAL,'noholes');
        BCOORDINATES=B{1};%get boundary coordinates
        %assemble boundary image
        BOUNDARY=zeros(size(SIGNAL));
        for i=1:length(BCOORDINATES)
            BOUNDARY(BCOORDINATES(i,1),BCOORDINATES(i,2))=1;
        end
        %add boundary to BORDER
        BORDER(BOUNDARY==1)=1;
        %shrink SIGNAL by new BOUNDARY
        SIGNAL(BOUNDARY==1)=0;
    end
end
%replace BORDER pixels with unfiltered data
for i=1:size(BORDER,1)
    for j=1:size(BORDER,2)
        if BORDER(i,j)==1
            MEDC(i,j)=C(i,j);
        end
    end
end
%replace matrix by filtered matrix
D=MEDC;
%restore SIGNAL matrix
SIGNAL=SIGNALCOPY;
%add bad pixels (C=0) to SIGNAL matrix
SIGNAL(C==0)=0;
%% Data range
%determine data range
DLIST=[];
for i=1:size(D,1)
    for j=1:size(D,2)
        if SIGNAL(i,j)==1 && D(i,j)>0
            DLIST=[DLIST;[i,j,D(i,j)]];
        end
    end
end
[B,IX]=sort(DLIST(:,3));
SORTDLIST=zeros(size(DLIST));
for i=1:length(IX)
    SORTDLIST(i,:)=DLIST(IX(i),:);
end
%% map data in [0 1] interval for color coding
MAPD=zeros(size(D));
%Normalize DATA matrix
for i=1:size(D,1)
    for j=1:size(D,2)
        if SIGNAL(i,j)>0
            if D(i,j)>highlimit
                MAPD(i,j)=1;
            elseif lowlimit<=D(i,j) && D(i,j)<=highlimit
                MAPD(i,j)=(D(i,j)-lowlimit)/(highlimit-lowlimit);
            end
        end
    end
end
%% Resize data matrix
%signals matrix
    RSIGNAL=imresize(SIGNAL,resizefactor,'bicubic');
    RSIGNAL(RSIGNAL<0.5)=0;RSIGNAL(RSIGNAL>=0.5)=1;
    %data matrix
    RMAPD=imresize(MAPD,resizefactor,'bicubic');
    %filter data matrix
    RMAPD(RSIGNAL<1)=0;%apply RESSIGNALSB matrix
    RMAPD(RMAPD>1)=1;    
    [IMAPD,mapMAPD]=gray2ind(RMAPD,256);
    RGBMAPD=ind2rgb(IMAPD,jet(256));
    %color bad pixels (C=0) in white
    for i=1:size(RGBMAPD,1)
        for j=1:size(RGBMAPD,2)
            if RSIGNAL(i,j)<1
                RGBMAPD(i,j,:)=[1,1,1];
            end
        end
    end
%background image
if exist('SHIFTB')==1
    BADJUST=imadjust(SHIFTB);
    [IB,mapB]=gray2ind(BADJUST,256);%convert background to indexed image
    RGBIMAGE=ind2rgb(IB,gray(256));
end
%% create figure and axis
if exist('SHIFTB')==1
    overlayfig=figure('Name',stackfile,'MenuBar','none','Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],'Color','w','Visible','off');
    overlayaxes=axes('Units','normalized','Position',[0 0 1 1],'Visible','on','Drawmode','fast');    
    
    heartfig=figure('Name',stackfile,'MenuBar','none','Units','normalized',...
    'Position',[0.25 0.25 0.5 0.5],'Color','w','Visible','off');
    heartaxes=axes('Units','normalized','Position',[0 0 1 1],'Visible','on','Drawmode','fast');
end

mapfig=figure('Name',file,'MenuBar','none','Units','normalized',...
'Position',[0.25 0.25 0.5 0.5],'Color','w','Visible','off');
mapaxes=axes('Units','normalized','Position',[0 0 1 1],'Visible','on','Drawmode','fast');


%draw images
axes(mapaxes)
image(RGBMAPD,'alphadata',RSIGNAL);

if exist('SHIFTB')==1
    %overlay image
    axes(overlayaxes)
    image(RGBIMAGE);hold on
    image(RGBMAPD,'alphadata',RSIGNAL*0.65);
    %heart image
    axes(heartaxes)
    image(RGBIMAGE)
end
set(mapaxes,'TickLength',[0 0])
set(mapfig,'visible','on');
set(mapaxes,'YDIR','reverse');
set(mapaxes,'DataAspectRatio',[1 1 1]);
set(mapaxes,'Visible','off');

if exist('SHIFTB')==1
    set(heartaxes,'TickLength',[0 0])
    set(heartfig,'visible','on');
    set(heartaxes,'YDIR','reverse');

    set(overlayaxes,'TickLength',[0 0])
    set(overlayfig,'visible','on');
    set(overlayaxes,'YDIR','reverse');
end

%draw scalebars at the bottom right corner of map image
if showscalebar==1
    scalebarlength=unitlength/pixelcalfactor;%length of scalebar in pixels (80x80 image)
    fractionscalebar=0.10;
    unitstring=[num2str(unitlength),' \mum'];
    scalebarxdata=resizefactor*(80-scalebarlength+1):resizefactor*fractionscalebar:resizefactor*80;
    scalebarydata=resizefactor*80*ones(1,length(scalebarxdata));
    xoffset=-10;
    yoffset=-20;
    axes(mapaxes)
    line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','k','LineWidth',2);
    line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','k','LineWidth',2);
    %add unit length to horizontal scalebar
    %text('Position',[resizefactor*70+xoffset+20,resizefactor*78+yoffset],'String',unitstring,'color','k','FontWeight','bold','FontSize',10)
end

%% add frame
%% select zoom area and plot new contours
%define width of frame
axes(mapaxes)

    %select corners of frame
    if setframeumwidth==0 || setframeumheight==0
        fprintf('select two corners of frame \n');
        for k=1:2
            [ax(k),ay(k)]=ginput(1);
            %rectangle('Position',[ax(k)-0.5,ay(k)-0.5,1,1],'LineWidth',2,'EdgeColor','r')
        end    
        FrameXLim=round([min(ax),max(ax)]);FrameYLim=round([min(ay),max(ay)]);
        framewidth=FrameXLim(2)-FrameXLim(1);frameheight=FrameYLim(2)-FrameYLim(1);
        frameumwidth=round(framewidth*pixelcalfactor);frameumheight=round(frameheight*pixelcalfactor);
        FRAMECENTER=[mean(FrameYLim),mean(FrameXLim)];
    else
        %select center of frame
        if exist('FRAMECENTER')==1 && isempty(FRAMECENTER)==0
            ax=FRAMECENTER(2);ay=FRAMECENTER(1);
        else
            fprintf('select center of frame \n');
            [ax,ay]=ginput(1);ax=round(ax/resizefactor);ay=round(ay/resizefactor);
            FRAMECENTER(1)=ay;FRAMECENTER(2)=ax;
        end
        framewidth=setframeumwidth/pixelcalfactor;frameheight=setframeumheight/pixelcalfactor;    
        FrameXLim=round([ax,ax+framewidth]-framewidth/2);
        FrameYLim=round([ay,ay+frameheight]-frameheight/2);
        
        
        %check if Frame limits are consistent with array limit
        if FrameXLim(1)<1,framewidth=framewidth-(1-FrameXLim(1));FrameXLim(1)=1;end
        if FrameXLim(2)>size(D,2),framewidth=framewidth-(size(D,2)-FrameXLim(2));FrameXLim(2)=size(D,2);end
        if FrameYLim(1)<1,frameheight=frameheight-(1-FrameYLim(1));FrameYLim(1)=1;end
        if FrameYLim(2)>size(D,1),frameheight=frameheight-(size(D,1)-FrameYLim(2));FrameYLim(2)=size(D,1);end        
        
        frameumwidth=(FrameXLim(2)-FrameXLim(1))*pixelcalfactor; frameumheight=(FrameYLim(2)-FrameYLim(1))*pixelcalfactor;
    end

    
%show frame
if showframe==1
    rectangle('Position',resizefactor*[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
end
if showcenter==1
    rectangle('Position',resizefactor*[FRAMECENTER(2),FRAMECENTER(1),1,1],'LineWidth',2,'EdgeColor',framecolor)
end
if exist('NFRAME')>0 && sframe>0
    asksaveframe=input('save frame (1=yes)?');
    if asksaveframe==1
        NFRAME=[NFRAME;[FrameXLim,FrameYLim]];
    end
end
%% add colorbars
numcolors=225;
division=stepsize;
tickdivision=numcolors/highlimit*division;
colormap(jet(numcolors));
if showcolorbar==1
    axes(mapaxes)
    cbar=colorbar('peer',mapaxes,'position',colorbarposition,'YTickMode','manual','Box','on','YColor','k','XColor','k','FontWeight','bold','YAxisLocation','right','Layer','top');
    %set number of tickmarks
    YTick=[1;numcolors];
    YTickLabel={num2str(lowlimit);num2str(highlimit)};
    set(cbar,'YTick',YTick,'YTickLabel',YTickLabel);
end
drawnow;
%% export images
%construct file names
if exist('imagepath')==0 || isempty(imagepath)==1 || length(imagepath)<3
    imagepath=datapath;
end
%imagepath=uigetdir(imagepath,'pick path to save images');imagepath=[imagepath,'\'];
%description=input('DESCRIPTION>','s');

filemap=[file,'-',description,'.tif'];

fileheart=[stackfile(1:end-4),'-heart-',description,'.tif'];
fileoverlay=[stackfile(1:end-4),'-overlay-',description,'.tif'];
filetrans=[stackfile(1:end-4),'-trans-',description,'.txt'];
filedata=[stackfile(1:end-4),'-data-',description,'.txt'];


% filemap=[stackfile(1:end-4),'-map-',num2str(lowlimit),'-',num2str(highlimit),'-',description,'.tif'];
% fileheart=[stackfile(1:end-4),'-heart-',num2str(lowlimit),'-',num2str(highlimit),'-',description,'.tif'];
% fileoverlay=[stackfile(1:end-4),'-overlay-',num2str(lowlimit),'-',num2str(highlimit),'-',description,'.tif'];

%get current image and save as movie frame
%MAP
imagefile=[imagepath,filemap];
%axes(mapaxes)
set(mapaxes,'DataAspectratioMode','auto');
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%crop image to remove black border
IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
%IMAGEC=IMAGE;
EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
%write image to file
imwrite(EXPORTIMAGE,imagefile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
set(mapaxes,'DataASpectratioMode','manual');

if exist('SHIFTB')
    %HEART
    imagefile=[imagepath,fileheart];
    axes(heartaxes)
    F=getframe(gcf);
    %convert movieframe back into image
    [IMAGE,imagemap]=frame2im(F);
    %crop image to remove black border
    IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
    %IMAGEC=IMAGE;
    EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
    %write image to file
    imwrite(EXPORTIMAGE,imagefile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);

    %OVERLAY
    imagefile=[imagepath,fileoverlay];
    axes(overlayaxes)
    F=getframe(gcf);
    %convert movieframe back into image
    [IMAGE,imagemap]=frame2im(F);
    %crop image to remove black border
    IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
    %IMAGEC=IMAGE;
    EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
    %write image to file
    imwrite(EXPORTIMAGE,imagefile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
end
%% show results
%get mean transient
if limits<5
    T=[];A=[];B=[];FRAMEij=[];LOW=[];HIGH=[];AMP=[];totalscans=300;U=[];KS=[];DURij=[];
    totalscans=2000;
    threshold=0.5;
    for i=FrameYLim(1):FrameYLim(2)
        for j=FrameXLim(1):FrameXLim(2)
            %save stimulus for frame center
            if SIGNAL(i,j)>0
                %get transient
                if exist('EPISODECA')==1
                    A=squeeze(EPISODECA(i,j,:));
                else
                    if exist('EPISODEFITR')==1 && askfit==1
                        A=squeeze(EPISODEFITR(i,j,:));
                    else
                         A=squeeze(EPISODER(i,j,:));
                    end 
                end
                baseline=mean(A(squeeze(BSCANS(i,j,:))));
                %add data to make sure transients are long enough
                A=[ones(round(totalscans/2),1).*mean(A);A;ones(round(totalscans/2),1).*mean(A)];
                if exist('EPISODESTIM','var')==1 && isempty(EPISODESTIM)==0
                    ST=[ones(round(totalscans/2),1).*mean(EPISODESTIM);EPISODESTIM;ones(round(totalscans/2),1).*mean(EPISODESTIM)];
                    if i==FRAMECENTER(1) && j==FRAMECENTER(2)
                        S=[];%save one stimulus for FRAMECENTER signal
                        S=ST((k-round(totalscans/2)+1:k+round(totalscans/2)));
                    end
                end
                U=round(totalscans/2)+UPSTROKERANGE;
                %determine upstroke and align
                [mval,midx]=min(A(U));midx=midx+U(1)-1;
                U2=[midx-10:U(end)]';
                [val,idx]=max(A(U2));idx=idx+U2(1)-1;%maximum
                NA=A-baseline;NA=NA/NA(idx);%normalize transient
                k=idx;while NA(k)>threshold;k=k-1;end;k=k+1;%scan of 92% maximum
                B=A(k-round(totalscans/2)+1:k+round(totalscans/2));
                T=[T,B];
                %save index to determine mean k to cut stimulus            
                KS=[KS;k];
                %get other measurements
                FRAMEij=[FRAMEij;[i,j]];
                if DUR(i,j)>0
                    DURij=[DURij;[DUR(i,j),i,j]];
                end
                if exist('CALOW')
                    LOW=[LOW;CALOW(i,j)];
                    HIGH=[HIGH;CAHIGH(i,j)];
                    AMP=[AMP;(CAHIGH(i,j)-CALOW(i,j))];
                else
                    LOW=[LOW;RLOW(i,j)];
                    HIGH=[HIGH;RHIGH(i,j)];
                    AMP=[AMP;(RHIGH(i,j)-RLOW(i,j))];
                end
            end
        end
    end
    %display results
    RESULTS=[[mean(LOW),std(LOW)];[mean(AMP),std(AMP)];[mean(DURij(:,1),1),std(DURij(:,1))]];
    p=[imagepath,filedata];
    %save(num2str(p),'RESULTS','-ascii');
    %modify stimulus to fit in data
    if exist('ST','var')==1
        MS=S/max(S);MS=MS*mean(AMP)/5;MS=MS-min(MS)+mean(LOW)+mean(AMP)/3;
    else
        MS=[];
    end
    TRANS=[([1:size(T,1)]'-round(totalscans)/2-1)/scanrate*1000,MS,mean(T,2)];
    p=[imagepath,filetrans];
    %save(num2str(p),'TRANS','-ascii');
    fprintf(['\n','mean baseline ratio: ',num2str(RESULTS(1,1)),' std: ',num2str(RESULTS(1,2)),'\n']);
    fprintf(['mean amplitude: ',num2str(RESULTS(2,1)),' std: ',num2str(RESULTS(2,2)),'\n']);
    fprintf(['mean 50pc duration: ',num2str(RESULTS(3,1)),' std: ',num2str(RESULTS(3,2)),'\n']);
    fprintf(['number of transients evaluated in selected area: ',num2str(size(T,2)),'\n']);
    
 elseif limits==5
    %get mean APD
    APDLIST=[];VMAXLIST=[];
    for i=FrameYLim(1):FrameYLim(2)
        for j=FrameXLim(1):FrameXLim(2)
            if SIGNAL(i,j)>0
                if APDMATRIX(i,j)>0
                    APDLIST=[APDLIST;[i,j,APDMATRIX(i,j)]];
                    VMAXLIST=[VMAXLIST;[i,j,VMAXMATRIX(i,j)]];
                end
            end
        end
    end
    fprintf(['\n','Mean APD80: ',num2str(mean(APDLIST(:,3))),' +-(SE) ',num2str(std(APDLIST(:,3))/sqrt(size(APDLIST,1))),' ms.\n']);
    fprintf(['\n','Mean Vmax: ',num2str(mean(VMAXLIST(:,3))),' +-(SE) ',num2str(std(VMAXLIST(:,3))/sqrt(size(VMAXLIST,1))),' s^-1.\n']);
end

%savedatabase:

if exist('DATABASE','var')==0
    DATABASE=[];id=1;
else
    id=size(DATABASE,2);
    size(DATABASE)
    fprintf(['previous heart, id:',num2str(id),' name:',num2str(DATABASE(id).name),'\n']);
    fprintf('continue? \n');
    %id=input('next id:');
    %hard code for efficiency (automatically continue previous id)
    id=1;
    %pause
end
%% general information
tracedownsample=1;%keep every nth sample
if id>length(DATABASE)
    generalinput=1;
else
    DATABASE(id)
    %generalinput=input('Edit general info? [RETURN]=NO, 1=YES>');
    %hard code to not edit general info for efficiency:
    generalinput=[];
    if isempty(generalinput)==1
        generalinput=0;
    end
end

if generalinput==1
    %DATABASE(id).age=input('age (hpf):');
    %DATABASE(id).name=stackfile;
    
    %nkx information
    %DATABASE(id).nkxtype=input('nkx mutant (homozygous ONLY) 0=TRUEWT 1=HET 25=nkx2.5 27=nkx2.7 2527=double:');
    
    %DATABASE(id).date=input('date (yymmdd):');
    %DATABASE(id).heart=input('heart number:');
    %DATABASE(id).rate=input('HR (bpm):');
    %DATABASE(id).scanrate=2000;
    %DATABASE(id).sinusloc=input('impulse initiation 0=atrium 1=av 2=ventricle 3=fstim:');
    %DATABASE(id).stim=input('stimulation rate 0:no stimulation (bpm):');
    
    %hard code inputs for efficiency:
    DATABASE(id).age=[];
    DATABASE(id).name=stackfile;
    DATABASE(id).date=[];
    DATABASE(id).heart=[];
    DATABASE(id).rate=[];
    DATABASE(id).scanrate=2000;
    DATABASE(id).sinusloc=[];
    DATABASE(id).stim=[];
    
    %rotate heart so that the atrium is on top. OFT is either left of
    %right. This defines the surface orientation of the heart
    %DATABASE(id).orientation=input('With A on top, OFT is LEFT(1) or RIGHT(2):');
    DATABASE(id).orientation=[];
    
    %stretch information
%   DATABASE(id).stretchregion=input('Area stretched (1) Atrium, (2) Ventricle:');
%   DATABASE(id).stretchtime=input('stretch time (minutes):');
%   DATABASE(id).recovertime=input('recovery time (minutes):');

    
    %DATABASE(id).comment=input('comment:','s');
    DATABASE(id).comment=[];
    
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
    if exist('VPLOTBIN','var')==1
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
if exist('VPLOTBIN', 'var') == 1
    % Continue only, if conduction data are available
fprintf(['current id:',num2str(id),' \n']);
DATABASE(id).name
fprintf('continue? \n');
%auto-enter for efficiency
[]
ESC

%region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
%hard code for efficiency, =DESCRIPTION
region=description

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
% if isfield(DATABASE,'ANGLEUNITVECTORXY')==1 && (isempty(DATABASE(id).ANGLEUNITVECTORXY)==1 || sum(DATABASE(id).ANGLEUNITVECTORXY==angleunitvector)<2)
%    DATABASE(id).ANGLEUNITVECTORXY=angleunitvector;%angle unit vector for specific region
%    fprintf('Unitvector changed.\n');pause
% end


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

%newname=input(['Name of database [',num2str(lastname),']>'],'s');
%auto-enter for efficiency
newname=lastname

if isempty(newname)==1
    filename=[lastname,'.mat'];
else
    filename=[newname,'.mat'];
end
lastname=filename(1:end-4);

%vinput=input(['WRITE TO FILENAME:',num2str(filename),' [ENTER]'],'s');
%auto-enter for efficiency
vinput=num2str(filename)
save(filename,'DATABASE');

end
%% SAVE APDs
if exist('APDMATRIX')==1
    if exist('vinput','var')==0
        %region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
        region=description
    end
    if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') ...
        || strcmp(region,'sa')==1 || strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
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
            if exist('NSTIM','var')==1 && isempty(NSTIM)==0 && i==FrameYLim(1) && j==FrameXLim(1)
                ST=[ADDSCANS;NSTIM;ADDSCANS];
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
if exist('ST','var')==1 && isfield(DATABASE,'STIM')==1
    if isempty(DATABASE(id).STIM)==1
        DATABASE(id).STIM=ST(round(mean(KLIST))-beforescans:round(mean(KLIST))+afterscans);
    end
end

%DATABASE(id).DATAscanrate=scanrate/tracedownsample;
tracedownsample=DATABASE(id).scanrate/DATABASE(id).DATAscanrate;

%DATABASE(id).(regionframedata)=downsample(T,tracedownsample);
DATABASE(id).(regionframedata)=T;
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
%auto-enter for efficiency
newname=lastname

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