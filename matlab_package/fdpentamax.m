%filter optical mapping data
function [FDATA]=fdpentamax(RDATA,lpcutoff,spatpix,scanrate)
%% determine if RDATA is ordered in colums or frames
s=size(RDATA);
FDATA=double(zeros(size(RDATA)));%allocate memory
%% filter coefficients
%[A42,B42]=butter(6,[36 48]/(scanrate/2),'stop');
[A28,B28]=butter(6,[25 30]/(scanrate/2),'stop');
%low pass filter
[A,B]=butter(6,lpcutoff/(scanrate/2),'low');% lp filter
%% FILER ALGORITHMS
if length(s)==3
    %TEMPORAL PROCESSLING
    hdl = waitbar(0,['IMAGE STACK TEMPORAL PROCESSING']);
    for j=1:size(RDATA,2)
        for i=1:size(RDATA,1)
            %extract time series
            PIXEL=double(squeeze(RDATA(i,j,:)));
            %mirror data
            MPIXEL=[flipud(PIXEL(1:200));PIXEL;flipud(PIXEL(end-200:end))];
            %FILTERS
            N_PIXEL=filtfilt(A28,B28,MPIXEL-mean(MPIXEL))+mean(MPIXEL);
            LP_PIXEL=filtfilt(A,B,N_PIXEL-mean(N_PIXEL))+mean(N_PIXEL);
            %median filter
            MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),5)+mean(LP_PIXEL);    
            FDATA(i,j,:)=MED_PIXEL(201:200+length(PIXEL));
        end
        waitbar(j/size(RDATA,2));
    end
    close(hdl);
    %SPATIAL PROCESSLING
    hdl = waitbar(0,['IMAGE STACK SPATIAL PROCESSING']);
    %spatial filter loop
    for frame=1:size(FDATA,3)
        WIEN=wiener2(FDATA(:,:,frame),[spatpix,spatpix]);
        MFILT=medfilt2(WIEN,[spatpix,spatpix]);
        FDATA(:,:,frame)=MFILT;
        waitbar(frame/size(FDATA,3));
    end
    close(hdl);
end