function [FM] = fdpolmatrix(M,spatw,spatm)
%fdpolmatrix: apply spatial filters wiener2 and medfilt2 to time matrices
%M:time matrix, spatpix:pixels, SIGNAL

%% filter algorithm
SPATWIEN=[spatw spatw];
SPATMED=[spatm spatm];

%WIENER FILTER
if spatw==0
    WIENM=M;
else
    WIENM=wiener2(M,SPATWIEN);
end

%MEDIAN FILTER
if spatm==0
    MEDM=WIENM;
else
    MEDM=medfilt2(WIENM,SPATMED);
end

%identify data boundary
SIGNAL=ones(size(M));
SIGNAL(M==0)=0;
%take off one boundary pixel
B1=bwperim(SIGNAL);SIGNAL1=SIGNAL;SIGNAL1(B1==1)=0;
B2=bwperim(SIGNAL1);SIGNAL2=SIGNAL1;SIGNAL2(B2==1)=0;
%take off boundary pixel from filtered matrix
if spatw==0 && spatm==0
    FM=M;
else
    FM=MEDM;FM(SIGNAL1==0)=0;
end
