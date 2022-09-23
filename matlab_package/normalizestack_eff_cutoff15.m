%clear RAWDATA FILTERDATA
%ONLY REQUIRED FOR 32 BIT OS!
%buildstack;%->pathfile
%filterstack;%->FILTERDATA
%bleachfit;%->BLEACHDATA
maxsn=1.5;
%% load data
if exist('BLEACHDATA')==0 || isempty(BLEACHDATA)==1
    [loadfile datapath]=uigetfile('-B.mat','load BLEACHDATA');
    loadpathfile=[datapath,loadfile];
    load(num2str(loadpathfile));
end
%% define search interval
if exist('PIXELCOORDINATES')==0
    %pick sample pixels
    %[PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,SCATTER,ASIGNALPIXELS,DATA)
    if exist('FBR')==0
        FBR=[];
    end
    [PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,[],[],FILTERDATA);
end
%% select upstroke range and baseline
if exist('EPISODETRACE','var')==0
for i=1:size(PIXELCOORDINATES,1)
    EPISODETRACE(:,i)=squeeze(BLEACHDATA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
end
end

if exist('UPSTROKERANGE')==0
    %plot TRACES
    figure('name','UPSTROKE RANGE')
    plot(EPISODETRACE);hold on
    fprintf('select range of upstrokes (2 clicks)\n');
    [ax,ay]=ginput(2);
    UPSTROKERANGE=[round(ax(1)):round(ax(2))]';
    plot(UPSTROKERANGE,EPISODETRACE(UPSTROKERANGE,:),'m');
    plot(1:length(EPISODETRACE),0,'k');
end
%% ask to repeat analysis if ASIGNALPIXELS available
if exist('ASIGNALPIXELS','var')==1
    askanalysis=input('USE EXISTING SIGNAL MATRIX [RETURN]');
    if isempty(askanalysis)==0
        ASIGNALPIXELS=[];
    end
end
if exist('ASIGNALPIXELS','var')==0 || isempty(ASIGNALPIXELS)==1
%% Exclude background from search
defaultcutoff=0.1;
if exist('FIRSTFRAME','var')==0
if exist(num2str(stackpath),'dir')==0
    [framefile framepath]=uigetfile('*.tif','Load image from high speed camera');
    FIRSTFRAME=imread([framepath,framefile],'tif');
else
    [framefile framepath]=uigetfile('*.tif','Load image from high speed camera or cancel',stackpath);
    if length(framefile)<2
        FIRSTFRAME=imread([stackpath,stackfile(1:end-4),'-FL.tif']);%load fluorescence image
    else
        FIRSTFRAME=imread([framepath,framefile],'tif');
    end
end
end

if exist('RGBFBR')==0
    FBR2=imadjust(mat2gray(FBR));
    [FBRIDX,map]=gray2ind(FBR2,256);
    RGBFBR=ind2rgb(FBRIDX,gray(256));
end

NEWFIRSTFRAME=FIRSTFRAME;%new first frame with replaced pixels

if exist('cutoff')==0
        cutoff=defaultcutoff;
end

NEWFIRSTFRAME=FIRSTFRAME;%new first frame with replaced pixels

h=figure;
repeat=1;

while repeat==1
    repeat=0;%do not repeat
    figure(h);imagesc(NEWFIRSTFRAME);
    
    lthreshold=cutoff*(max(FIRSTFRAME(:))-min(FIRSTFRAME(:)));%threshold for the detection of Background intensities
    DARKPIXELS=ones(size(BLEACHDATA,1),size(BLEACHDATA,2));%ARRAY of DARKPIXELS (binary image)
    NEWFIRSTFRAME=FIRSTFRAME;

    for i=1:size(FIRSTFRAME,1)
        for j=1:size(FIRSTFRAME,2)
            if FIRSTFRAME(i,j)<lthreshold
                DARKPIXELS(i,j)=0;%mark background pixel in binary image
                NEWFIRSTFRAME(i,j)=0;
            end
        end
    end

    %SHOW BACKGROUND CUTOFF
    figure(h);
    IFRAME=mat2gray(NEWFIRSTFRAME);
    [XFRAME,mapXFRAME]=gray2ind(IFRAME,256);%convert to index
    RGBFRAME=ind2rgb(XFRAME,jet(256));%convert to RGB
    ALPHA=alphamatrix(NEWFIRSTFRAME,0.2,0.7);
    image(RGBFBR);hold on;image(RGBFRAME,'AlphaData',ALPHA);hold off

    %Repeat if image not satisying
    %newcutoff=input(['Enter background cutoff [%] of range [',num2str(cutoff*100),']>']);
    newcutoff=15 %defaults to 15
    %if isempty(newcutoff)==0
        cutoff=newcutoff/100;
    %   repeat=1;    
    %end
end
%% find signals
SN=zeros(size(DARKPIXELS));
MAXIDX=zeros(size(SN));
MAXVAL=zeros(size(SN));
%% search loop
hdl = waitbar(0,['DETERMINE SIGNAL-TO-NOISE RATIOS']);
for i=1:size(BLEACHDATA,1)
    for j=1:size(BLEACHDATA,2)  
    %%signal to noise ratios
    if DARKPIXELS(i,j)==1            
        %extract data
        PIXEL=-double(squeeze(BLEACHDATA(i,j,:)));
        %sort values
        UPSTROKEAMP=[UPSTROKERANGE,PIXEL(UPSTROKERANGE)];
        SORTUPSTROKEAMP=sortrows(UPSTROKEAMP,2);
        %average the largest amplitudes for 5 ms
        averagen=round(0.005*scanrate)+1;
        MAXVAL(i,j)=mean(SORTUPSTROKEAMP(end-averagen:end,2));
        %get good index for maximum
        MAXIDX(i,j)=mean(SORTUPSTROKEAMP(end-averagen:end,1));
        %get noise data
        if exist('BASELINESCANS','var')==1
            BASELINE=squeeze(BASELINESCANS(i,j,:));
        end
        NOISEDATA=PIXEL(BASELINE);
        %sort baseline scans and measure amplitude
        SORTNOISEDATA=sort(NOISEDATA);   
        noiseamplitude=mean(SORTNOISEDATA(end-averagen:end))-mean(SORTNOISEDATA(1:averagen));
        %estimate S/N
        SN(i,j)=(MAXVAL(i,j)-mean(PIXEL(BASELINE)))/noiseamplitude;
    else
       SN(i,j)=0;
    end
    %% new pixel
    end
    waitbar(i/size(BLEACHDATA,1));
end
close(hdl)
%% BUILD BINARY IMAGE BASED ON SN
ASIGNALPIXELS=zeros(size(BLEACHDATA,1),size(BLEACHDATA,2));
ASIGNALPIXELS(SN>=maxsn)=1;%signal-to-noise ratio
I1=ASIGNALPIXELS;
I2=I1;
%I2=imfill(I1,'holes');
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
end
%% display binary image of SIGNALPIXELS (ASIGNALPIXELS)
bw=figure('Name','ASIGNALPIXELS','MenuBar','none','Units','normalized','Position',[0 0.3 0.5 0.5],'Color','w','Visible','on');
bwaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);
tr=figure('Name','TRACE','MenuBar','none','Units','normalized','Position',[0.5 0.3 0.5 0.5],'Color','w','Visible','on');

axes(bwaxes);
imshow(ASIGNALPIXELS,'InitialMagnification','fit');
% button=0; commented out for efficiency (auto "ESC")

while button~=27
    figure(bw);
    [ax,ay,button]=ginput(1);
    if button~=27
        i=round(ay);j=round(ax); 
        %plot PIXEL
        BASELINE=squeeze(BASELINESCANS(i,j,:));
        PIXEL=-double(squeeze(BLEACHDATA(i,j,:))); 
        %need to measure maximum amplitude if selected pixel outside
        %DARKPIXELS
        if DARKPIXELS(i,j)==0
            UPSTROKEAMP=[UPSTROKERANGE,PIXEL(UPSTROKERANGE)];
            SORTUPSTROKEAMP=sortrows(UPSTROKEAMP,2);
            MAXVAL(i,j)=mean(SORTUPSTROKEAMP(end-averagen:end,2));
        end
        
        maxpixel=MAXVAL(i,j)-mean(PIXEL(BASELINE));
        PIXEL=PIXEL-mean(PIXEL(BASELINE));
        PIXEL=PIXEL/maxpixel;
        
        %mark 50% activation threshold
        figure(tr);
        plot(PIXEL,'k');
        hold on;
            plot(1:length(PIXEL),0.5,'r');
            plot(BASELINE,PIXEL(BASELINE),'b.');
            plot(1:length(PIXEL),0,'b');
            plot(UPSTROKERANGE,PIXEL(UPSTROKERANGE),'r');
        hold off
        drawnow
        %add position to the list of pixels to be transferred
        if ASIGNALPIXELS(i,j)==0
            ASIGNALPIXELS(i,j)=1;    
            %mark pixel in binary image
            figure(bw)
            rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
            drawnow;
        else
            figure(bw)
            ASIGNALPIXELS(i,j)=0;
            rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','k'); 
            drawnow;
        end
    end
end
%% show result
h=figure('Name','SIGNALS');
imshow(im2bw(ASIGNALPIXELS,0.5),'InitialMagnification','fit');
%% ASSEMBLE NORMDATA
clear NORMDATA DATA RAWDATA
NORMDATA=zeros(size(BLEACHDATA));
hdl = waitbar(0,['normalizing DATA']);
for i=1:size(ASIGNALPIXELS,1)
    for j=1:size(ASIGNALPIXELS,2)
        if ASIGNALPIXELS(i,j)==1
            PIXEL=-double(squeeze(BLEACHDATA(i,j,:)));              
            BASELINE=squeeze(BASELINESCANS(i,j,:));
            maxpixel=MAXVAL(i,j)-mean(PIXEL(BASELINE));
            PIXEL=PIXEL-mean(PIXEL(BASELINE));
            PIXEL=PIXEL/maxpixel;
            %save PIXEL in NORMDATA
            NORMDATA(i,j,:)=PIXEL;%write PIXEL back into DATA array
        end
    end
    waitbar(i/size(ASIGNALPIXELS,1))
end
close(hdl)
%% normalize stimulus
stimbaseline=min(EPISODESTIM(:));%define stimulus baseline
OSTIM=EPISODESTIM-stimbaseline;%remove stimulus offset
NSTIM=OSTIM/max(OSTIM(:));
%% save NORMDATA array
if exist('datapath')==0
    [newfile datapath]=uiputfile('*.mat','DATA FOLDER');
end

[fpath, fname, fext] = fileparts(stackfile);
normfile=[stackfile(1:end-length(fext)),'-N.mat'];

normpathfile=[datapath,normfile];
fprintf(['SAVING DATA in \n'])
normpathfile
if exist('STIM')==0
    STIM=[];
end

if exist('SEARCHMAX')==0
    SEARCHMAX=[1,size(NORMDATA,3)];
end

save(normpathfile,'NSTIM','NORMDATA','EPISODEFRAME','EPISODEX',...
    'BASELINESCANS','ASIGNALPIXELS','BASELINE','UPSTROKERANGE',...
    'PIXELCOORDINATES','SN','MAXIDX','MAXVAL','DARKPIXELS',...
    'scanrate','datapath','normfile','stackfile',...
    'firstframenumber','lastframenumber','FR','RGBFBR','lpcutoff','spatpix');
