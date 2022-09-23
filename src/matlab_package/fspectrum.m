%create fluorescence excitation spectrum with monochromator
%CCD A(:,1), wav A(:,2), STIM=A(:,3)
%% load data

if exist('file')==0
    [file path]=uigetfile('*.txt','LOAD FLUORESCENCE DATA');
    pathfile=[path file];
    CCD=load(pathfile);

    [wfile wpath]=uigetfile('*.txt','LOAD WAVELENGTH DATA');
    wpathfile=[wpath wfile];
    WAV=load(wpathfile);
end
A=[CCD,WAV];
%% scan parameter
fprintf('DATA=[CCD,DAC]\n')
scanrate=input(['scanrate [SCANS/S]:']);
steptime=input('time per step [MS]:');
stepscans=steptime*1.0e-3*scanrate;
%% filter data
%low pass filter
FLDATA=A(:,1);%fluorescence data
WDATA=A(:,2);%wavelength data

[FLA,FLB]=butter(6,100/(scanrate/2),'low');% lp filter fluorescence
[WA,WB]=butter(6,100/(scanrate/2),'low');% lp filter BC output

%median filter
M_FLDATA=medfilt1(FLDATA-mean(FLDATA),5)+mean(FLDATA);
M_WDATA=medfilt1(WDATA-mean(WDATA),5)+mean(WDATA);
%M_WDATA=WDATA;

%low-pass filter
FM_FLDATA=filtfilt(FLA,FLB,M_FLDATA-mean(M_FLDATA))+mean(M_FLDATA);
FM_WDATA=filtfilt(WA,WB,M_WDATA-mean(M_WDATA))+mean(M_WDATA);

%% detect scan boundaries
SPEC=FM_FLDATA;WAV=FM_WDATA;

STARTIDX=[];ENDIDX=[];

minwav=min(WAV);maxwav=max(WAV);

idx=11;
while idx<length(WAV)-10
    idx=idx+1;
    if WAV(idx-1)>minwav+0.7*(maxwav-minwav) && WAV(idx)<minwav+0.7*(maxwav-minwav) && WAV(idx+10)<minwav+0.3*(maxwav-minwav);
        STARTIDX=[STARTIDX;idx+2];
    end
end

%% prepare data for display and select spectrum
%correct offsets
if exist('A(:,3)')==1
    STIM=A(:,3);
else
    STIM=zeros(size(SPEC));
end
M=[SPEC,WAV,STIM];
for j=1:size(M,2)
    M(:,j)=M(:,j)-min(M(:,j));
    M(:,j)=M(:,j)/max(M(:,j));
end
    
%plot data
figure
plot(M)

%%select spectrum
ax=1;
%[ax,ay,b]=ginput(1);

%find data selected
if ax<=STARTIDX(1)%selected first spectrum
        XDATA=1:STARTIDX(1)-1;
else
    for i=1:length(STARTIDX)-1
        if ax>STARTIDX(i) && ax<=STARTIDX(i+1)
            XDATA=(STARTIDX(i):STARTIDX(i+1)-1)';
        end
    end
end


%transform voltage to wavlength data
W=WAV(XDATA)/10+300;
S=SPEC(XDATA);
%% set center scan points
dstep=round(stepscans/2-5);%number of scans to average left and right of center wavelength

i=length(W)-round(stepscans/2-1);%scans in W

CENTERSCAN=[];
while i>2*dstep
    CENTERSCAN=[CENTERSCAN;i];%scans in center of W
    i=i-stepscans;
end

%re-order CENTER scans wavelendth 
CENTERSCAN=flipud(CENTERSCAN);
CENTERTIMESCAN=CENTERSCAN-length(W)+XDATA(end);



%calculate wavelength and intensity for each CENTERSCAN

%plot spectrum with calculated ranges
figure
plot(S,'k');
hold on

LAMBDA=zeros(length(CENTERSCAN),1);
INTENSITY=zeros(size(LAMBDA));
for i=1:length(CENTERSCAN)
    LAMBDA(i)=mean(W(CENTERSCAN(i)-dstep:CENTERSCAN(i)+dstep));
    INTENSITY(i)=mean(S(CENTERSCAN(i)-dstep:CENTERSCAN(i)+dstep));
    plot(CENTERSCAN(i),S(CENTERSCAN(i)),'r+');
    plot(CENTERSCAN(i)-dstep:CENTERSCAN(i)+dstep,INTENSITY(i),'m-');
end

%fit linear function to wavelength data
YDATA=LAMBDA;
XDATA=(1:length(LAMBDA))';
P=polyfit(XDATA,YDATA,1);%LINEAR MODEL
FITLAMBDA=polyval(P,XDATA);



%fit spectrum with cubic spline
LAMBDA=FITLAMBDA;%use fitted wavelength
splinesteps=1000;
lambdastepsize=(LAMBDA(end)-LAMBDA(1))/splinesteps;
scanstepsize=(CENTERTIMESCAN(end)-CENTERTIMESCAN(1))/splinesteps;

SPLINE=csapi(LAMBDA,INTENSITY);
SPLINELAMBDA=(LAMBDA(1):lambdastepsize:LAMBDA(end))';%new wavelength axis
SPLINEINTENSITY=fnval(SPLINE,SPLINELAMBDA);
SPLINESCAN=(CENTERTIMESCAN(1):scanstepsize:CENTERTIMESCAN(end))';%time axis


%plot spectrum
figure
plot(LAMBDA,INTENSITY,'r+');
hold on
plot(SPLINELAMBDA,SPLINEINTENSITY,'k')


%% tramsform wavelength axis into 'scans after stimulus' axis

    if exist('A(:,3)')==1
        STIM=A(:,3);
        %normalize stimulus
        NSTIM=STIM-mean(STIM(round(ax)-10:round(ax)+10),1);
        NSTIM=NSTIM/max(NSTIM);

        %find index of stimulus
        stimidx=round(ax);
        threshold=0.5;
        while NSTIM(stimidx)<0.5 && stimidx>1,stimidx=stimidx-1;end%falling edge
        while NSTIM(stimidx)>0.5 && stimidx>1,stimidx=stimidx-1;end%rising edge
    else
        stimidx=0;
    end

    %transform CENTERTIMESCAN axis
    CENTERTIMESCAN=CENTERTIMESCAN-stimidx;
    SPLINESCAN=SPLINESCAN-stimidx;
%% save results
filepoints=['POINTS-',file(1:end-4),'.txt'];
filespline=['SPLINE-',file(1:end-4),'.txt'];

% POINTS=[LAMBDA,CENTERTIMESCAN,INTENSITY];
% SPLINE=[SPLINELAMBDA,SPLINESCAN,SPLINEINTENSITY];

POINTS=[LAMBDA,INTENSITY];
SPLINE=[SPLINELAMBDA,SPLINEINTENSITY];


save([path,filepoints],'POINTS','-ascii','-tabs');
save([path,filespline],'SPLINE','-ascii','-tabs');