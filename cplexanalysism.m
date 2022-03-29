%process cardioplex data

scanrate=2000;

if exist('RAWDATA')==0
    [file path]=uigetfile('*.txt','CARDIOPLEX ASCII DATA FILE');
    RAWDATA=load([path,file],'-ascii');
    RAWDATA=double(RAWDATA);
    [stimfile stimpath]=uigetfile('*.txt','STIMULUS OR CANCEL');
    if stimfile~=0
        RAWSTIM=load([stimpath,stimfile],'-ascii');
        %filter stimulus
        [ASTIM,BSTIM]=butter(6,125/(scanrate/2),'low');% lp filter
        %low-pass filter
        LP_STIM=filtfilt(ASTIM,BSTIM,RAWSTIM-mean(RAWSTIM))+mean(RAWSTIM);
        %median filter
        STIM=medfilt1(LP_STIM-mean(LP_STIM),40)+mean(LP_STIM);
    else
        RAWSTIM=[];STIM=[];
    end
end
%% filter fluorescence data
FILTERDATA=double(zeros(size(RAWDATA)));%allocate memory

%NOTCH FILTER1
[A500,B500]=butter(6,[480 520]/(scanrate/2),'stop');

%NOTCH FILTER2
[A150,B150]=butter(6,[100 200]/(scanrate/2),'stop');

%NOTCH FILTER3
[A030,B030]=butter(5,[25 35]/(scanrate/2),'stop');

%low pass filter
[A,B]=butter(6,100/(scanrate/2),'low');% lp filter

for i=1:size(RAWDATA,2)
    PIXEL=double(RAWDATA(:,i));
    %mirror data
    MPIXEL=[flipud(PIXEL(1:200));PIXEL;flipud(PIXEL(end-200:end))];
    %notch filters
    N_PIXEL=filtfilt(A500,B500,MPIXEL-mean(MPIXEL))+mean(MPIXEL);
    NN_PIXEL=filtfilt(A150,B150,N_PIXEL-mean(N_PIXEL))+mean(N_PIXEL);
    NNN_PIXEL=NN_PIXEL;%SKIPT 3rd notch filter
    %NNN_PIXEL=filtfilt(A030,B030,NN_PIXEL-mean(NN_PIXEL))+mean(NN_PIXEL);
    %LP filter   
    LP_PIXEL=filtfilt(A,B,NNN_PIXEL-mean(NNN_PIXEL))+mean(NNN_PIXEL);
    %median filter
    MED_PIXEL=medfilt1(LP_PIXEL-mean(LP_PIXEL),40)+mean(LP_PIXEL);
    %wavelet filter
        %calculate default parameters
        [C,L]=wavedec(MED_PIXEL,3,'db10');
        %default parameter for de-noising
        [thr,sorh,keepapp]=ddencmp('den','wv',MED_PIXEL);
        %signal reconstruction
        W_PIXEL=wdencmp('gbl',C,L,'db10',3,thr,sorh,keepapp);
        %save processed data
   FILTERDATA(:,i)=W_PIXEL(201:200+length(PIXEL));
end
%% isolate action potential with stimulus
trace=figure('Name','CLICK LEFT AND RIGHT OF PEAK');
plot(-FILTERDATA(:,1));

[AX,AY]=ginput(2);
XDATA=round([AX(1):AX(2)]);
YDATA=-FILTERDATA(XDATA,:);

if isempty(STIM)==0
    YSTIM=STIM(XDATA);
end    

%% correct photobleaching
bfig=figure('Name','select 2 baseline intervals (4 clicks, left to right)');
plot(YDATA(:,1));
[ax,ay]=ginput(4);
ax=round(ax);
BSCANS=[ax(1):ax(2),ax(3):ax(4)]';
NEWYDATA=[];
XDATA=[1:size(YDATA,1)]';
for i=1:size(YDATA,2)
    P=polyfit(BSCANS,YDATA(BSCANS,i),1);%LINEAR MODEL
    NEWYDATA(:,i)=YDATA(:,i)-((polyval(P,1:length(XDATA))')-mean(YDATA(ax(1):ax(2),i),1));
    %correct offset
    NEWYDATA(:,i)=NEWYDATA(:,i)-mean(NEWYDATA(ax(1):ax(2),i));
    %normalize and align
    [val,idx]=max(NEWYDATA(:,i));
    NEWYDATA(:,i)=NEWYDATA(:,i)/mean(NEWYDATA(idx-10:idx+10,i));
end
%process stimulus
if exist('YSTIM')==1
    YSTIM=YSTIM-mean(YSTIM);
    %normalize and align
    NEWSTIM=YSTIM/max(YSTIM);
end
CUTDATA=NEWYDATA;
CUTSTIM=NEWSTIM;
DCUTDATA=diff(CUTDATA);
nfig=figure('Name','normalized action potential and stimulus');
plot(CUTDATA,'r');hold on;plot(CUTSTIM,'k');
%% calculate APD
thresholdleft=.5;
thresholdright=.5;
for i=1:size(CUTDATA,2)
%% process column signal i
N=CUTDATA(:,i);
[val,maxidx]=max(N);
%find index to 10% on both sides
left=maxidx;
while N(left)>thresholdleft && left>1;left=left-1;end;
    %better estimate by interpolation
    LEFTX=[left-2:left+3];LEFTY=N(LEFTX);
    LEFTXi=[left-2:0.01:left+3];LEFTYi=interp1(LEFTX,LEFTY,LEFTXi);
    k=1;while LEFTYi(k)<thresholdleft;k=k+1;end;left2=LEFTXi(k);

right=maxidx;
while N(right)>thresholdright && right<length(N);right=right+1;end;
    %better estimate by interpolation
    RIGHTX=[right-10:right+10];RIGHTY=N(RIGHTX);
    RIGHTXi=[right-10:0.01:right+10];RIGHTYi=interp1(RIGHTX,RIGHTY,RIGHTXi);
    k=1;while RIGHTYi(k)>thresholdright;k=k+1;end;right2=RIGHTXi(k);

%show APD calculation
h=figure('Name',['AP ',num2str(i)]);
plot(N);hold on
xplot=[left2:0.1:right2];
plot(xplot,thresholdright,'r');
%calculate APD in ms
APD(i)=(right2-left2)/scanrate*1.0e3;
%% loop
end
%% collect results
%AP and DAP arrays
DAP=[zeros(1,size(CUTDATA,2));DCUTDATA]*scanrate;
if exist('CUTSTIM')==1
    AP=[CUTSTIM,CUTDATA];
else
    AP=CUTDATA;
end

%Vmax in s^-1
k=figure('Name','CLICK LEFT AND RIGHT OF PEAK');plot(DAP);
[dx,dy]=ginput(2);

file
APD
D=max(DAP(round(dx(1)):round(dx(2)),:))

% save data
save('AP.txt','AP','-ascii');
save('DAP.txt','DAP','-ascii');