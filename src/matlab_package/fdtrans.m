function [FDATA]=fdtrans(RDATA,lpcutoff,medscans,scanrate)
% filter transient data
% data needs to be organized in columns
S=size(RDATA);
FDATA=double(zeros(size(RDATA)));%allocate memory
%filter coefficients
[A,B]=butter(6,lpcutoff/(scanrate/2),'low');% lp filter
%% FILER ALGORITHMS
if length(S)==2 && S(2)<S(1) 
    %TEMPORAL PROCESSLING
    for i=1:size(RDATA,2)
        %extract time series
        PIXEL=RDATA(:,i);
        %mirror data
        M_PIXEL=[flipud(PIXEL(1:200));PIXEL;flipud(PIXEL(end-200:end))];
        %LP filter   
        LP_PIXEL=filtfilt(A,B,M_PIXEL-mean(M_PIXEL))+mean(M_PIXEL);
        %median filter
        MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),medscans)+mean(LP_PIXEL);    
        %save trace
        FDATA(:,i)=MED_PIXEL(201:200+length(PIXEL));
    end
end