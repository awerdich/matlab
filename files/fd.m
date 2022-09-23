%filter optical mapping data
function [FDATA]=fd(RDATA,lpcutoff,spatpix,median,scanrate)
%% determine if RDATA is ordered in colums or frames
TDATA=double(zeros(size(RDATA)));%allocate memory to store temporally filtered data
FDATA=zeros(size(TDATA));%allocate memory to store spatially filtered data;
%% filter coefficients BAND PASS
[A500,B500]=butter(6,[450 540]/(scanrate/2),'stop');
[A120,B120]=butter(6,[90 140]/(scanrate/2),'stop');
%LOW PASS
[A,B]=butter(6,lpcutoff/(scanrate/2),'low');% lp filter
%% FILTER ALGORITHMS
tstart=tic;%start timer
%TEMPORAL PROCESSLING
hdl = waitbar(0,['IMAGE STACK TEMPORAL PROCESSING']);
for i=1:size(RDATA,1)
    parfor j=1:size(RDATA,2)
        %extract time series
        PIXEL=double(squeeze(RDATA(i,j,:)));
        %mirror data
        MPIXEL=[flipud(PIXEL(1:200));PIXEL;flipud(PIXEL(end-200:end))];
        %notch filters
        N_PIXEL=filtfilt(A500,B500,MPIXEL-mean(MPIXEL))+mean(MPIXEL);
        NN_PIXEL=filtfilt(A120,B120,N_PIXEL-mean(N_PIXEL))+mean(N_PIXEL);
        %LP filter   
        LP_PIXEL=filtfilt(A,B,NN_PIXEL-mean(NN_PIXEL))+mean(NN_PIXEL);
        %median filter
        MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),median)+mean(LP_PIXEL);    
        %save trace
        TDATA(i,j,:)=MED_PIXEL(201:200+length(PIXEL));
    end
    waitbar(i/size(RDATA,1));
end
close(hdl);

%SPATIAL PROCESSLING OPERATING ON TDATA
%spatial filter loop
parfor frame=1:size(TDATA,3)
    F=squeeze(TDATA(:,:,frame));
    WIEN=wiener2(F,[spatpix,spatpix]);
    MED=medfilt2(WIEN,[spatpix,spatpix]);
    FDATA(:,:,frame)=MED;
end
elapsedtime=toc(tstart);fprintf(['Elapsed time: ',num2str(elapsedtime),' s\n']);