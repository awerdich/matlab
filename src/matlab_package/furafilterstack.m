%% filter parameter
showfilter=1;
lpcutoff=100;%cutoff frequency of lowpass filter
spatpix=3;%number of pixels to include in spatial filtering in each direction
scanrate=250;%interpolated scanrate [ratios/s]
DATAST(1).reverselam=0;%set to 1 to switch wavelengths
camrate=DATAST(1).scanrate;
stackpath=DATAST(1).PATH;
stackfile=DATAST(1).FILE;
fratio=DATAST(1).fratio;
PLOTLOCATION=[40,40];
stimchannel=2;%BNC channel for stimulus data
%% sort wavelengths in stack
[RW1,RW2]=furasortwavelengths(DATAST);
%reverse wavelengths if flag set in DATAST
if DATAST(1).reverselam==1
    RW3=RW1;RW1=RW2;RW2=RW3;clear RW3;
end
%write fluorescence image
FBR=imadjust(mat2gray(RW1(:,:,20)));
[FBRIDX,map]=gray2ind(FBR,256);
DATAST(1).FBR=FBR;
DATAST(1).RGBFBR=ind2rgb(FBRIDX,gray(256));
%% load calibration data
%load calibration data to filter constants
% if exist('calfile')==0
%     [calfile calpath]=uigetfile('*.mat','LOAD CALIBRATION RESULTS OR CANCEL');
% end
% if length(calfile)>1
%     load([calpath,calfile],'KD','F0','RMIN','RMAX');
% end
%% interpolate data to 250 ratios/s
if (camrate/fratio)<scanrate
    %interpolate
    hdl = waitbar(0,['INTERPOLATING IMAGE STACK']);
    RT=((1:size(RW1,3))-1)*fratio/camrate;%original time base
    T=RT(1):1/scanrate:RT(end);%new time base
    W1=zeros(size(RW1,1),size(RW1,2),length(T));%first wavelength
    W2=zeros(size(RW2,1),size(RW2,2),length(T));%second wavelength
    for i=1:size(RW1,1)
        for j=1:size(RW2,2)
            RPIXELW1=squeeze(RW1(i,j,:));
            W1(i,j,:)=interp1(RT,RPIXELW1,T,'pchip');
            RPIXELW2=squeeze(RW2(i,j,:));
            W2(i,j,:)=interp1(RT,RPIXELW2,T,'pchip');
        end
        waitbar(i/size(W1,1));
    end
close(hdl)
else
    W1=RW1;
    W2=RW2;
end
%% load stimulus
RAWSTIM=DATAST(1).BNCDATA(:,stimchannel);
%correct offset and normalize
OSTIM=RAWSTIM-mean(RAWSTIM);
NSTIM=OSTIM/max(OSTIM);
%interpolate or downsample to fit ratio data
%downsample data
STIM =downsample(NSTIM,round(camrate/scanrate));
STIM=STIM(1:size(W1,3)); 

% Normalize stimulus
STIM_O = STIM - mean(STIM);
STIM_N = STIM_O/max(STIM_O);
STIM = STIM_N;
%% filter coefficients
%low pass filter
[A,B]=butter(6,lpcutoff/(scanrate/2),'low');% lp filter
%% filter one pixel
if showfilter==1
    %extract time series
    PIXEL=double(squeeze(W2(PLOTLOCATION(1),PLOTLOCATION(2),:)));
    %mirror data
    MPIXEL=[flipud(PIXEL(1:200));PIXEL;flipud(PIXEL(end-200:end))];
    %lowpass
    LP_PIXEL=filtfilt(A,B,MPIXEL-mean(MPIXEL))+mean(MPIXEL);
    %median filter
    MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),20)+mean(LP_PIXEL);
    %wavelet filter
        %calculate default parameters
        [C,L]=wavedec(MED_PIXEL,3,'db10');
        %default parameter for de-noising
        [thr,sorh,keepapp]=ddencmp('den','wv',MED_PIXEL);
        %signal reconstruction
        W_PIXEL=wdencmp('gbl',C,L,'db10',3,thr,sorh,keepapp);
        T=W_PIXEL;%signal after temporal filtering
    
    %calculate spectra
    [P_PIXEL,f_PIXEL]=pwelch(MPIXEL-mean(MPIXEL),256,[],[],scanrate);
    [P_WPIXEL,f_WPIXEL]=pwelch(W_PIXEL-mean(W_PIXEL),256,[],[],scanrate);

    %plot spectra and filtered data
    figure
    plot(f_PIXEL,P_PIXEL,'k');hold on
    plot(f_WPIXEL,P_WPIXEL,'r')
    figure
    plot(MPIXEL(201:200+length(PIXEL)),'k');hold on
    plot(LP_PIXEL(201:200+length(PIXEL)),'b');hold on
    plot(MED_PIXEL(201:200+length(PIXEL)),'r');
    plot(W_PIXEL(201:200+length(PIXEL)),'m');
end
%% filter whole stack in time
FW1=zeros(size(W1));%first wavelength
FW2=zeros(size(W2));%second wavelength
hdl = waitbar(0,['IMAGE STACK TEMPORAL PROCESSING']);
for j=1:size(W1,2)
    parfor i=1:size(W1,1)
    %extract time series
    PIXELW1=double(squeeze(W1(i,j,:)));
    PIXELW2=double(squeeze(W2(i,j,:)));
    %mirror data
    MPIXELW1=[flipud(PIXELW1(1:200));PIXELW1;flipud(PIXELW1(end-200:end))];
    MPIXELW2=[flipud(PIXELW2(1:200));PIXELW2;flipud(PIXELW2(end-200:end))];
    %LP filter   
    LP_PIXELW1=filtfilt(A,B,MPIXELW1-mean(MPIXELW1))+mean(MPIXELW1);
    LP_PIXELW2=filtfilt(A,B,MPIXELW2-mean(MPIXELW2))+mean(MPIXELW2);
    %median filter
    MED_PIXELW1=medfilt1(LP_PIXELW1-mean(LP_PIXELW1),20)+mean(LP_PIXELW1);
    MED_PIXELW2=medfilt1(LP_PIXELW2-mean(LP_PIXELW2),20)+mean(LP_PIXELW2);
    %wavelet filter
    %calculate default parameters
    [CW1,LW1]=wavedec(MED_PIXELW1,3,'db10');
    [CW2,LW2]=wavedec(MED_PIXELW2,3,'db10');
    %default parameters for de-noising
    [thrW1,sorhW1,keepappW1]=ddencmp('den','wv',MED_PIXELW1);
    [thrW2,sorhW2,keepappW2]=ddencmp('den','wv',MED_PIXELW2);
    %signal reconstruction
    W_PIXELW1=wdencmp('gbl',CW1,LW1,'db10',3,thrW1,sorhW1,keepappW1);
    W_PIXELW2=wdencmp('gbl',CW2,LW2,'db10',3,thrW2,sorhW2,keepappW2);
    %save processed data
    FW1(i,j,:)=W_PIXELW1(201:200+length(PIXELW1));
    FW2(i,j,:)=W_PIXELW2(201:200+length(PIXELW2));
    end
    waitbar(j/size(W1,2));
end
close(hdl);
%% calculate ratios
lam1=DATAST(1).lam1;
lam2=DATAST(1).lam2;
if exist('KD')==1
    
    FR=zeros(size(FW1));%RATIOS
    FCA=zeros(size(FW1));%CALIBRATED RATIOS
    BACKGROUND1=zeros(size(FW1));
    BACKGROUND2=zeros(size(FW1));
    
    hdl = waitbar(0,['APPLY CALIBRATION DATA']);
    for i=1:size(FW1,1)
        for j=1:size(FW1,2)
            PIXEL1=squeeze(FW1(i,j,:));
            PIXEL2=squeeze(FW2(i,j,:));         
            %ratio calculation
            if lam1<lam2
                R=(PIXEL1-BACKGROUND1(i,j))./(PIXEL2-BACKGROUND2(i,j));
            else
                R=(PIXEL2-BACKGROUND2(i,j))./(PIXEL1-BACKGROUND1(i,j));
            end
            FCA(i,j,:)=KD(i,j)*F0(i,j)*(R-RMIN(i,j))./(RMAX(i,j)-R);
            FR(i,j,:)=R;
        end
        waitbar(i/size(W1,1));
    end
    close(hdl);
else
    %calculate ratios without background information
    hdl = waitbar(0,['CALCULATE RATIOS']);
    for i=1:size(FW1,1)
        for j=1:size(FW1,2)
            PIXEL1=squeeze(FW1(i,j,:));
            PIXEL2=squeeze(FW2(i,j,:));         
            if lam1<lam2
                R=PIXEL1./PIXEL2;
            else
                R=PIXEL2./PIXEL1;
            end
            FR(i,j,:)=R;
        end
        waitbar(i/size(W1,1));
    end
    close(hdl);
end

% apply spatial filtering
if spatpix>0
    hdl = waitbar(0,['IMAGE STACK SPATIAL PROCESSING']);
    for i=1:size(FR,3)
        FRAME=wiener2(FR(:,:,i),[spatpix,spatpix]);
        FR(:,:,i)=FRAME;
        if exist('FCA')==1
            FRAME=wiener2(FCA(:,:,i),[spatpix,spatpix]);
            FCA(:,:,i)=FRAME;
        end
        waitbar(i/size(FR,3));
    end
    close(hdl)
end
clear FRAME
furableachfit