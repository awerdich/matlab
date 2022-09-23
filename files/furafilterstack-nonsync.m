%% read image stack
if exist('scanrate')==0
    scanrate=input('scanrate>');
end

%check if filename and path available was loaded or re-load if not
%available
if exist('pathstackfile')==0
    if exist('stackpath')==1 && exist('stackfile')==1
    pathstackfile=[stackpath,stackfile];
    end
    [stackfile stackpath]=uigetfile('*.tif','MULTI-FRAME TIF FILE');
    pathstackfile=[stackpath,stackfile];
end
if exist('RAWDATA')==0
    info=imfinfo(num2str(pathstackfile));
    nframes=size(info,2);
    RAWDATA=uint16(zeros(info(1).Width,info(1).Height,nframes)); %pre-allocate 3D array 
    %read images into DATA array
    hdl = waitbar(0,['reading image stack',num2str(stackfile)]);
    for frame=1:nframes
        [RAWDATA(:,:,frame)]=imread(pathstackfile,frame);
        waitbar(frame/nframes);
    end
    close(hdl)
end
%% read wavelength and ready signals
if exist('WAV')==0
    [wavfile wavpath]=uigetfile('*.txt','WAVELENGTH DATA FOR ALL FRAMES');
    PREWAV=load([wavpath,wavfile],'-ascii');
    [readyfile readypath]=uigetfile('*.txt','READY DATA FOR ALL FRAMES OR CANCEL');
    if readyfile~=0
        PREREADY=load([readypath,readyfile],'-ascii');
    else
        PREREADY=[];
    end
    %get frame numbers in tif file
    if exist('firstframenumber')==0 || exist('lastframenumber')==0
        firstframenumber=input('first frame number>');
        lastframenumber=input('last frame number>');
    end
    %cut wavelength data to corresponding fluorescence data 
    WAV=PREWAV(firstframenumber:lastframenumber);
    if isempty(PREREADY)==0
        READY=PREREADY(firstframenumber:lastframenumber);
    else
        READY=[];
    end
end
%% load stimulus
[stimfile stimpath]=uigetfile('*.txt','STIMULUS OR CANCEL');
if stimfile~=0
    ALLSTIM=load([stimpath stimfile],'-ascii');
    %cut stimulus
    RAWSTIM=ALLSTIM(firstframenumber:lastframenumber);   
    %filter stimulus
    [ASTIM,BSTIM]=butter(6,125/(scanrate/2),'low');% lp filter
    %low-pass filter
    LP_STIM=filtfilt(ASTIM,BSTIM,RAWSTIM-mean(RAWSTIM))+mean(RAWSTIM);
    %median filter
    MED_STIM=medfilt1(LP_STIM-mean(LP_STIM),40)+mean(LP_STIM);
    %correct offset and normalize
    STIM=MED_STIM-mean(MED_STIM(1:200));
    STIM=STIM/max(STIM);
end
%% get wavelengths and corresponding scan intervals
[WAVELENGTHS,SCANINTERVALS]=findwav(WAV,READY);
showseparation;
%% generate time basis for each wavlength
XDATA=[];
for k=1:size(SCANINTERVALS,1)
    for m=1:length(WAVELENGTHS)
        XDATA(k,m)=mean(SCANINTERVALS(k,2*m-1:2*m));
    end
end
%% determine average scanrate and set filter coefficients
SCAN=[];
for lambda=1:length(WAVELENGTHS)
    SCAN(lambda)=1/(mean(diff(XDATA(:,lambda)))/scanrate);
end
MEANSCAN=mean(SCAN);
%low pass filter coefficients
[A,B]=butter(6,80/(MEANSCAN/2),'low');% lp filter
%% separate wavelengths for each pixel time series and filter time series
%determine maximum number of frames
FDATAR=zeros(size(RAWDATA,1),size(RAWDATA,2),size(XDATA,1));
FDATA1=zeros(size(FDATAR));%rawdata at first wavelength
FDATA2=zeros(size(FDATAR));%rawdata at second wavelength
hdl = waitbar(0,'filtering fluorescence ratios');
for i=1:size(RAWDATA,1)
    for j=1:size(RAWDATA,2)
%% calculate mean fluorescence data for PIXEL
        PIXEL=double(squeeze(RAWDATA(i,j,:)));
        F=zeros(size(XDATA,1),length(WAVELENGTHS));
        R=zeros(size(XDATA,1),1);
        SPIXEL=zeros(size(XDATA,1),3);
        for k=1:size(XDATA,1)
            for lambda=1:length(WAVELENGTHS)
                I=[SCANINTERVALS(k,2*lambda-1),SCANINTERVALS(k,2*lambda)];
                X=(I(1):I(2));
                F(k,lambda)=mean(PIXEL(X));
            end
        end
        R=F(:,1)./F(:,2);
        SPIXEL=[F(:,1),F(:,2),R];
        %filter fluorescence data
        FPIXEL=zeros(size(SPIXEL));
        if min(SPIXEL(:,1))>0
            for k=1:size(SPIXEL,2)
                %mirror data
                MR=[flipud(SPIXEL(1:10,k));SPIXEL(:,k);flipud(SPIXEL(end-10:end,k))];
                %low-pass filter
                LP=filtfilt(A,B,MR);
                %median filter
                MED=medfilt1(LP,4);
                %wavelet filter
                %calculate default parameters
                [C,L]=wavedec(MED,3,'db10');
                %default parameter for de-noising
                [thr,sorh,keepapp]=ddencmp('den','wv',MED);
                %signal reconstruction
                WLT=wdencmp('gbl',C,L,'db10',3,thr,sorh,keepapp);
                %remove excess data
                FPIXEL(:,k)=WLT(11:10+size(SPIXEL,1));
            end
        end
%% loop
    FDATA1(i,j,:)=FPIXEL(:,1);
    FDATA2(i,j,:)=FPIXEL(:,2);
    FDATAR(i,j,:)=FPIXEL(:,3);
    end
    waitbar(i/size(RAWDATA,1))
end
close(hdl)
%% Spatial filter
hdl = waitbar(0,['spatial filter']);
% spatial filter loop
DATA=FDATAR;
FILTERDATA=DATA;
spatinterval=[8 8];
for frame=1:size(DATA,3)
    WIEN=wiener2(DATA(:,:,frame),spatinterval);   
    MFILT=medfilt2(WIEN,spatinterval);
    FILTERDATA(:,:,frame)=MFILT;
    waitbar(frame/size(DATA,3));
end
close(hdl);
%% save data
filterfile=[stackfile(1:end-4),'-F.mat'];
filterpath=stackpath;%save mat file in tif folder
filterpathfile=[filterpath,filterfile];
fprintf(['saving FILTERDATA in \n'])
filterpathfile
pause;
if exist('RAWSTIM')==1
    save(filterpathfile,'RAWDATA','RAWSTIM','STIM','FILTERDATA','scanrate','SCAN','SCANINTERVALS','WAV','PREWAV','READY','PREREADY','WAVELENGTHS','SCANINTERVALS','XDATA','stackfile','stackpath','filterfile','filterpath','firstframenumber','lastframenumber');
else
    save(filterpathfile,'RAWDATA','FILTERDATA','scanrate','SCAN','SCANINTERVALS','WAV','PREWAV','READY','PREREADY','WAVELENGTHS','SCANINTERVALS','XDATA','stackfile','stackpath','filterfile','filterpath','firstframenumber','lastframenumber');
end