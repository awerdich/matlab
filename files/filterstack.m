%process stack from single file
%uses pathstackfile returned by buildstack.m
%read raw data into 3D array
%pre-allocate memory
showfilter=1;
PLOTLOCATION=[40,40];
lpcutoff=120;%cutoff frequency of lowpass filter CAREFUL, AFFECTS CVs!!!!
spatpix=4;%number of pixels in spatial filter 2 for good S/N 4 for bad S/N3
medianfilter=40;%median filter
invertdata=1;%invert data if peaks go up
%% get image stack name and location
if exist('DATAST','var')==1
    stackpath=DATAST(1).PATH;
    stackfile=DATAST(1).FILE;
    pathstackfile=[stackpath,stackfile];
    RAWDATA=DATAST(1).RAWDATA(:,:,10:end); %cut the first 10 frames
    DARKFRAME=DATAST(1).DARKFRAME;
    scanrate=DATAST(1).scanrate;
    firstframenumber=1;
    lastframenumber=size(RAWDATA,3);
else
    if exist('stackpath')==1 && exist('stackfile')==1
        pathstackfile=[stackpath,stackfile];
    end
    [stackfile stackpath]=uigetfile('*.tif','MULTI-FRAME TIF FILE');
    pathstackfile=[stackpath,stackfile];
    firstframenumber=input('first frame number>');
    lastframenumber=input('last frame number>');
    scanrate=input('scanrate>');
end
%% load stimulus
RAWSTIM=[];
if exist('DATAST','var')==1
    RAWSTIM=DATAST(1).BNCDATA(:,1);
else
    [stimfile stimpath]=uigetfile('*.txt','STIMULUS TEXT FILE OR CANCEL',stackpath);
    if stimfile~=0
        %import data
        IMSTRUCT=importdata([stimpath,stimfile],',',1);
        ALLSTIM=double(IMSTRUCT.data);
        %cut stimulus
        RAWSTIM=ALLSTIM(firstframenumber:lastframenumber);   
    end
end
if isempty(RAWSTIM)==0
    %filter stimulus
    [ASTIM,BSTIM]=butter(6,125/(scanrate/2),'low');% lp filter
    %low-pass filter
    LP_STIM=filtfilt(ASTIM,BSTIM,RAWSTIM-mean(RAWSTIM))+mean(RAWSTIM);
    %median filter
    MED_STIM=medfilt1(LP_STIM-mean(LP_STIM),40)+mean(LP_STIM);
    %correct offset and normalize
    STIM=MED_STIM-mean(MED_STIM(1:200));
    STIM=STIM/max(STIM);
else
    STIM=[];
end
%% load image stack
%check if filename and path available was loaded or re-load if not
%available
if exist('RAWDATA')==0
    info=imfinfo(num2str(pathstackfile));
    nframes=size(info,1);
    RAWDATA=double(zeros(info(1).Height,info(1).Width,nframes)); %pre-allocate 3D array 
    %read images into DATA array
    hdl = waitbar(0,['reading ',num2str(stackfile)]);
    for frame=1:nframes
        [RAWDATA(:,:,frame)]=imread(pathstackfile,frame);
        waitbar(frame/nframes);
    end
    close(hdl)
end
%% Normalize RAWDATA array by DARKFRAME
%correct for gain differences between recordings
%pixel gains can vary between recordings; 
%DARKFRAME is the average of 10 frames before shutter opening
%DARKFRAME is an OFFSET and NOT GAIN!!!!
if exist('DARKFRAME','var')==0
    [darkfile darkpath]=uigetfile('*.tif','OPEN DARK FRAME',stackpath);
    DARKFRAME=double(imread([darkpath,darkfile]));
end
if exist('FR','var')==0
    for i=1:size(RAWDATA,3)
        RAWDATA(:,:,i)=RAWDATA(:,:,i)-DARKFRAME;
    end
    FR=imadjust(mat2gray(RAWDATA(:,:,1)));
    imwrite(FR,[stackpath,stackfile(1:end-4),'-FL.tif']);    
end
%% invert data if needed
if invertdata==1
    RDATA=invertstack(RAWDATA);
else
    RDATA=RAWDATA;
end
%% FILTER FUNCTION
[FILTERDATA]=fd(RDATA,lpcutoff,spatpix,medianfilter,scanrate);
%% plot result
if showfilter==1
    %re-scale filtered data
    
    RPIXEL=squeeze(RDATA(PLOTLOCATION(1),PLOTLOCATION(2),:));
    FPIXEL=squeeze(FILTERDATA(PLOTLOCATION(1),PLOTLOCATION(2),:));

    %plot rawdata and filtered data
    figure
    plot(RPIXEL-mean(RPIXEL),'k');hold on
    plot(FPIXEL-mean(FPIXEL),'r');
        
    %spectra
    [P_RPIXEL,F_RPIXEL]=pwelch(RPIXEL-mean(RPIXEL),256,[],[],scanrate);
    [P_FPIXEL,F_FPIXEL]=pwelch(FPIXEL-mean(FPIXEL),256,[],[],scanrate);

    %plot spectra 
    figure;
    plot(F_RPIXEL,P_RPIXEL,'k');hold on   
    plot(F_FPIXEL,P_FPIXEL,'r')
    set(gca,'YLim',[0 10]);
    
end
%%  clear variables
clear DATA FDATA PIXEL FLIPPIXEL FFLIPPIXEL FFFLIPPIXEL