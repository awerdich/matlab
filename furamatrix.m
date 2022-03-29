%DETERMINE FURA MATRICES FOR CA RELEASE AND BASELINE CA CONCENTRATIONS
threshold=0.5;minsn=3;%minimum SNR
DATA=DATAST(1).EPISODER;
%% define background pixels
%check if ASIGNALPIXELS exist
startsearch=1;%start signal search if ASIGNALPIXELS not available
if exist('ASIGNALPIXELS')==1
   searchsignals=input('USE PREVIOUS SIGNAL ARRAY [RETURN]');
   if isempty(searchsignals)==1
       startsearch=0;%do not seach for signals, only normalize
   end
else
    searchsignals=0;
end
if startsearch==1
%% define search interval
if exist('PIXELCOORDINATES')==0
    if exist('stackpath')==0
        stackpath=uigetdir('FOLDER FOR IMAGE STACK'); 
        stackpath=[stackpath,'\'];
    end
    FBR=imread([stackpath,stackfile(1:end-4),'-FL.tif']);%load fluorescence image
    [PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,[],[],FR);
end
%% select upstroke range and baseline
if exist('UPSTROKERANGE')==0
    figure('name','range of baseline+upstrokes ')
    SPIXELS=[squeeze(DATA(PIXELCOORDINATES(1,2),PIXELCOORDINATES(1,1),:)),squeeze(DATA(PIXELCOORDINATES(2,2),PIXELCOORDINATES(2,1),:))];
    %correct offset and normalize
    for i=1:size(SPIXELS,2)
        SPIXELS(:,i)=SPIXELS(:,i)-mean(SPIXELS(:,i));
        SPIXELS(:,i)=SPIXELS(:,i)/max(SPIXELS(:,i));
    end
    plot(SPIXELS(:,1),'k');hold on
    plot(SPIXELS(:,2),'r');
    fprintf('select range baseline + upstrokes\n');
    [ax,ay]=ginput(2);
    UPSTROKERANGE=[round(ax(1)):round(ax(2))]';
    plot(UPSTROKERANGE,SPIXELS(UPSTROKERANGE,:),'m');
end
%% Exclude background from search
defaultcutoff=0.1;
%load bfrightfield image of heart for overlay
RGBFRAME=DATAST(1).RGBFBR;
FRAME=FBR;
if exist('cutoff')==0
        cutoff=defaultcutoff;
end


h=figure;
repeat=1;

while repeat==1
    repeat=0;%do not repeat
    figure(h);image(RGBFRAME);
    
    lthreshold=cutoff*(max(FRAME(:))-min(FRAME(:)));%threshold for the detection of Background intensities
    DARKPIXELS=ones(size(FRAME));%ARRAY of DARKPIXELS (binary image)

    NEWFRAME=zeros(size(FRAME));
    DARKPIXELS=zeros(size(FRAME));
    NEWFRAME(FRAME>=lthreshold)=1;
    DARKPIXELS(FRAME>=lthreshold)=1;
    
    %SHOW BACKGROUND CUTOFF
    figure(h);
    [INEWFRAME,mapINEWFRAME]=gray2ind(NEWFRAME,256);%convert to index
    RGBNEWFRAME=ind2rgb(INEWFRAME,jet(256));%convert to RGB
    ALPHA=alphamatrix(NEWFRAME,0.2,0.7);
    image(RGBFRAME);hold on;image(RGBNEWFRAME,'AlphaData',ALPHA);hold off

    %Repeat if image not satisying
    newcutoff=input(['Enter background cutoff [%] of range [',num2str(cutoff*100),']>']);

    if isempty(newcutoff)==0
        cutoff=newcutoff/100;
        repeat=1;    
    end
end

%% search loop
apixel=4;%number of pixels to average for amplitude measurement
baselinetime=175/2;%[ms] time of baseline before upstroke
baselinescans=floor(baselinetime*1.0e-3*scanrate)+1;%number of scans to average for measuring diastolic concentrations, baseline
maximumscans=10;%number of scans to average for measuring systolic concentrations
SN=zeros(size(DARKPIXELS));%signal-to-noise ratios
BSCANS=zeros(size(SN,1),size(SN,2),baselinescans);%baseline scans
MSCANS=zeros(size(SN,1),size(SN,2),maximumscans);%maximum scans
hdl = waitbar(0,['IDENTIFYING SINGNALS']);
for i=1:size(DATA,1)
    parfor j=1:size(DATA,2)  
%% Extract pixel       
        PIXEL=double(squeeze(DATA(i,j,:)));
        %measure amplitudes in UPSTROKERANGE
        UPSTROKEAMP=[];SORTUPSTROKEAMP=[];
        for k=UPSTROKERANGE(1)+2:UPSTROKERANGE(end);UPSTROKEAMP=[UPSTROKEAMP;[k-2,mean(PIXEL(k-round(apixel/2):k+round(apixel/2)-1))]];end
        SORTUPSTROKEAMP=sortrows(UPSTROKEAMP,2);
        %peak
        %delete early baseline scans to make sure maximum is measured
        %towards the end of the upstrokerange
        ASORTUPSTROKEAMP=SORTUPSTROKEAMP;
        D=[];
        for k=1:size(ASORTUPSTROKEAMP,1)
            if ASORTUPSTROKEAMP(k,1)<UPSTROKERANGE(1)+baselinescans
                D=[D;k];%mark index for deletion
            end
        end
        ASORTUPSTROKEAMP(D,:)=[];
        MAXIMUMSCANS=ASORTUPSTROKEAMP(end-maximumscans+1:end,1);
        MSCANS(i,j,:)=MAXIMUMSCANS;
        
        %baseline
        %delete sorted scans that appear after the maximum
        %to restrict baselinescans to interval before peak
        BSORTUPSTROKEAMP=SORTUPSTROKEAMP;
        D=[];
        for k=1:size(BSORTUPSTROKEAMP,1)
            if BSORTUPSTROKEAMP(k,1)>min(MAXIMUMSCANS)
                D=[D;k];%mark index for deletion
            end
        end
        BSORTUPSTROKEAMP(D,:)=[];
        
        %take find those scans that have the smallest amplitude before the
        %maximum
        B=BSORTUPSTROKEAMP(1:baselinescans,1);
        BALL=[min(B):max(B)]';
        BASELINESCANS=BALL(1:baselinescans);
        
        %take exaclty baselinescans after the first MINAMPBASELINESCAN
        BSCANS(i,j,:)=BASELINESCANS;
        
        %noise
        %define noisescans as the standard deviation of the baseline
        %maximum
        %noise=std(PIXEL(BASELINESCANS));
        noise=max(PIXEL(BASELINESCANS))-min(PIXEL(BASELINESCANS));
        
        %signal
        %calculate mean signal amplitude
        amp=mean(PIXEL(MAXIMUMSCANS))-mean(PIXEL(BASELINESCANS));
        %exclude signal if maximum was found before the first baselinescan
        if amp<0;amp=0;end
        %sn ratio
        SN(i,j)=amp/noise;
        
%% new pixel
    end
    waitbar(i/size(DATA,1));
end
close(hdl)
end
%% BUILD BINARY IMAGE BASED ON SN
ASIGNALPIXELS=zeros(size(DATA,1),size(DATA,2));
ASIGNALPIXELS(SN>=minsn)=1;%signal-to-noise ratio
ASIGNALPIXELS(DARKPIXELS<1)=0;

I1=ASIGNALPIXELS;
I2=imfill(I1,'holes');
% define objects in image
minpixel=100;%minimum number of pixels that make up an object
[LABEL,num]=bwlabel(I2,8);
%cound number of pixels in each object and delete all objects in I8 that
%are smaller than minpixel
DELOBJ=[];%[o-number,pixels] number of objects removed from image
for k=1:num
    %get pixel coordinates for object i
    [I,J] = find(LABEL==k);
    %delete object if it is too small
    if length(I)<minpixel
        for m=1:length(I)
            I2(I(m),J(m))=0;
        end
        DELOBJ=[DELOBJ;[k,length(I)]];%number of deleted object
    end
end
%re-count objects
[NEWLABEL,newnum]=bwlabel(I2,8);
%re-assign ASIGNALPIXELS
ASIGNALPIXELS=I2;
%% measure fura matrices
%data for ratios and concentrations
% apixel=4;%number of pixels to average for amplitude measurement
% baselinescans=20;%number of scans to average for measuring diastolic concentrations, baseline
% maximumscans=5;%number of scans to average for measuring systolic concentrations
% SN=zeros(size(DARKPIXELS));%signal-to-noise ratios
% BSCANS=zeros(size(SN,1),size(SN,2),baselinescans);%baseline scans
% MSCANS=zeros(size(SN,1),size(SN,2),maximumscans);%maximum scans
% NSCAN1=zeros(size(SN));%first noise scans
% NSCAN2=zeros(size(SN));%last noise scan

RLOW=zeros(size(ASIGNALPIXELS));%minimum R value
RHIGH=zeros(size(ASIGNALPIXELS));%maximum R value
LEFTIDX=zeros(size(ASIGNALPIXELS));%left index at threshold level
RIGHTIDX=zeros(size(ASIGNALPIXELS));%right index at threshold level
DUR=zeros(size(ASIGNALPIXELS));%transient duration (ms)

if exist('RCA')==1
CALOW=zeros(size(ASIGNALPIXELS));%minimum CA concentration
CAHIGH=zeros(size(ASIGNALPIXELS));%maximum CA concentration
end

hdl = waitbar(0,['CALCULATING FURA MATRICES']);
for i=1:size(ASIGNALPIXELS,1)
    for j=1:size(ASIGNALPIXELS,2)
        if ASIGNALPIXELS(i,j)>0
%% CALCULATE FURA MATRICES            
            %extract transient at location (i,j)
            PIXELR=squeeze(DATA(i,j,:));
            B=squeeze(BSCANS(i,j,:));
            M=squeeze(MSCANS(i,j,:));
            
            %allow measurements only on good signals
            if SN(i,j)>minsn
                 
                %define maximum and minumim and normalize transient
                if exist('CA','var')==1
                    PIXELCA=squeeze(CA(i,j,:));
                    BLIST=[B,PIXELCA(B)];%list of baseline amplitudes
                    SORTBLIST=sortrows(BLIST,2);%sort list to get smallest values
                    clow=mean(SORTBLIST(1:floor(baselinescans/2)+1,2));
                    chigh=mean(PIXELCA(M));
                    N=(PIXELCA-clow)/(chigh-clow);
                else
                    %take the baseline scans with the smallest amplitudes as
                    %diastolic level
                    BLIST=[B,PIXELR(B)];%list of baseline amplitudes
                    SORTBLIST=sortrows(BLIST,2);%sort list to get smallest values
                    rlow=mean(SORTBLIST(1:floor(baselinescans/2)+1,2));
                    rhigh=mean(PIXELR(M));
                    N=(PIXELR-rlow)/(rhigh-rlow);
                end
                %measure transient duration
                [upidx,downidx,NTRANS]=apd(N,threshold,UPSTROKERANGE,B,scanrate);
                
                %save data if duration measurement was successful
                if upidx>0 && downidx>0
                    RLOW(i,j)=rlow;
                    RHIGH(i,j)=rhigh;
                    if exist('CA')==1
                        PIXELCA=squeeze(CA(i,j,:));
                        CALOW(i,j)=clow;
                        CAHIGH(i,j)=chigh;
                    end
                    LEFTIDX(i,j)=upidx;
                    RIGHTIDX(i,j)=downidx;
                    DUR(i,j)=(downidx-upidx)/scanrate*1.0e3;
                end 
            end
            %% Next pixel
        end
    end
    waitbar(i/size(ASIGNALPIXELS,1));
end
close(hdl);
%% display binary image of SIGNALPIXELS (ASIGNALPIXELS)
bw=figure('Name','ASIGNALPIXELS','MenuBar','none','Units','pixels','Position',[10 250 500 500],'Color','w','Visible','on');
bwaxes=axes('Position',[0 0 1 1],'Visible','on','YDir','reverse','TickLength',[0 0]);
tr=figure('Name','TRACE','MenuBar','none','Units','pixels','Position',[510 250 500 500],'Color','w','Visible','on');
traxes=axes;
axes(bwaxes);
imshow(ASIGNALPIXELS,'InitialMagnification','fit');
button=0;
T=([1:size(DATA,3)]'-1)/scanrate;
while button~=27
    axes(bwaxes);
    [ax,ay,button]=ginput(1);
    if button~=27
        i=round(ay);j=round(ax); 
        %plot PIXEL
        PIXEL=double(squeeze(DATA(i,j,:))); 
        axes(traxes);plot(T,PIXEL,'k');
        B=squeeze(BSCANS(i,j,:));
        hold on
            plot(T(UPSTROKERANGE),PIXEL(UPSTROKERANGE),'m')
            plot(T(B),PIXEL(B),'g+')
            plot(T,RLOW(i,j),'b.');
            plot(T,RHIGH(i,j),'r.');
            ax=gca;
            %set(ax.XLabel,'String','Time (s)');
            %set(ax.YLabel,'String','Fluorescence ratio F340/F380');
        hold off
        drawnow
    end
end
%% show result
h=figure('Name','BASELINE');
imshow(im2bw(ASIGNALPIXELS,0.5),'InitialMagnification','fit');%process MOVETOSIGNAL list
%% save results
if exist('datapath')==0 || isempty(datapath)==1
    datapath=uigetdir(path,'FOLDER FOR FILTERED DATA'); 
    datapath=[datapath,'\'];
end
EPISODER=DATAST(1).EPISODER;
SCATTER=[];

[fpath, fname, fext] = fileparts(stackfile);
f2file=[stackfile(1:end-length(fext)),'-M.mat'];
durthreshold=threshold;
if exist('CA')==1
save([datapath,f2file],'UPSTROKERANGE','PIXELCOORDINATES','DARKPIXELS',...
'scanrate','datapath','f2file','stackfile',...
'ASIGNALPIXELS','SN','BSCANS','MSCANS','RLOW','RHIGH','CALOW','CAHIGH','LEFTIDX','RIGHTIDX','DUR','durthreshold');
else
 save([datapath,f2file],'UPSTROKERANGE','PIXELCOORDINATES','DARKPIXELS',...
'scanrate','datapath','f2file','stackfile',...
'ASIGNALPIXELS','SN','BSCANS','MSCANS','RLOW','RHIGH','LEFTIDX','RIGHTIDX','DUR','durthreshold',...
'EPISODER','SCATTER');
end
SCATTER=ones(size(ASIGNALPIXELS));