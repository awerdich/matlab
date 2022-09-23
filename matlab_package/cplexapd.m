%process cardioplex data
%% User parameter
scanrate=2000;
threshold=0.2;%amplitude threshold of normalized transient
totalscans=3000;%[ms] transient length
addscans=totalscans;
extremelength=40;% [ms] maximum and minimum length
escans=round(extremelength*1.0e-3*scanrate)+1;
beforemax=450;% [ms] transient length before upstroke
bscans=round(beforemax*1.0e-3*scanrate)+1;
clear RATE
calrate=0;%set to 1 if you want to estimate stimlation rates
%120 Hz notch filter
[A120,B120]=butter(8,[100 140]/(scanrate/2),'stop');
golaywindow=31;%smoothing filter window sent to Khaled on 10/6/2010 for cplexapd
%% open text file
if exist('textpath')==0 || length(textpath)<2
    textpath=[];
    textfile=[];
end
%get name and path
repeatloop=1;RAWDATA=[];i=1;
while repeatloop==1
    [textfile,textpath]=uigetfile('*.txt',['open CardioPlex text file. Last:',num2str(textfile)],textpath);%get filename
    if length(textfile)>1
        %import data
        IMSTRUCT=importdata([textpath,textfile],',',1);
        RAWDATA(i).DATA=IMSTRUCT.data;
        RAWDATA(i).FILE=textfile;
        RAWDATA(i).PATH=textpath;
        lastpath=textpath;i=i+1;
    else
        repeatloop=0;
        textpath=lastpath;
    end
end
stimcol=input('stimulus data in each file (0=no stimulus, 1=first, 2=last)>');
RATE=[];

if 0<stimcol && stimcol<=2
    %go through data and extract stimulus and frequency from each entry
    hdl = waitbar(0,['CALCULATING STIMULATION FREQUENCIES']);
    R=zeros(length(RAWDATA),1);
    NUMCOL=[];%count the number of colums in each file
    for i=1:length(RAWDATA)
        if stimcol==1
            RAWDATA(i).STIMULUS=RAWDATA(i).DATA(:,1);
            RAWDATA(i).DATA(:,1)=[];
        elseif stimcol==2
            RAWDATA(i).STIMULUS=RAWDATA(i).DATA(:,end);
            RAWDATA(i).DATA(:,end)=[];
        else
            RAWDATA(i).STIMULUS=[];
        end
        RAWDATA(i).FDATA=zeros(size(RAWDATA(i).DATA));        
        for j=1:size(RAWDATA(i).DATA,2);
            T=RAWDATA(i).DATA(:,j);
            FT=filtfilt(A120,B120,T-mean(T))+mean(T);
            RAWDATA(i).FDATA(:,j)=FT;
        end
        
    if calrate>0    
        %calcuate rate
        D=RAWDATA(i).STIMULUS;
        %set freqency vector [Hz]
        F=[0.01:0.001:5];
        %define spectrum object
        Hs=spectrum.periodogram;
        %options for spectrum object
        Hopts=msspectrumopts(Hs,(D-mean(D)));
        %set custom options
        Hopts.NormalizedFrequency=false;
        Hopts.Fs=scanrate;
        Hopts.SpectrumType='twosided';
        Hopts.FreqPoints='User Defined';
        Hopts.FrequencyVector=F;
        %calculate mean-squared spectrum
        %peaks reflect power in the signal at a given frequency band
        Hmss=msspectrum(Hs,(D-mean(D)),Hopts);
        %find maximum and normalize spectrum
        SPEC=Hmss.Data;
        [val,idx]=max(SPEC);NSPEC=SPEC/val;
        maxf=F(idx)*60;%heart rate in bpm
        RAWDATA(i).rate=maxf;
        RATE(i,1)=maxf;
    else
        RATE(i,1)=0;
    end
    waitbar(i/length(RAWDATA));
    %count columns in each file
    NUMCOL=[NUMCOL;size(RAWDATA(i).FDATA,2)];
    end
    close(hdl);
end 
fprintf('Stimulation frequencies:\n');
RATE
%select file i
%ALLAPD(i,j)=APD(j);
%ALLVMAX(i,j)=max(D1)*scanrate;


ALLAPD=zeros(length(RAWDATA),max(NUMCOL));
ALLVMAX=zeros(length(RAWDATA),max(NUMCOL));

DELDATA=[];%entries in APDATA to be deleted
for i=1:length(RAWDATA)
    repeat=1;
    while repeat==1
    %% get transients from file i
    TRANSD=-RAWDATA(i).FDATA;%invert voltage data
    if isfield(RAWDATA,'STIMULUS')==1
        STD=RAWDATA(i).STIMULUS;
        ST=STD-mean(STD);ST=ST/max(ST)/2;
        STIM=RAWDATA(i).STIMULUS;
        NSTIM=STIM-min(STIM);NSTIM=NSTIM/max(NSTIM);
        ASTIM=[zeros(addscans,size(STIM,2));STIM;zeros(addscans,size(STIM,2))];
    end
    %% define baseline and upstroke range
    fig=figure('NumberTitle','off');
    set(fig,'Name',['BASELINE UPSTROKE RANGE FILE: ',num2str(RAWDATA(i).FILE)]);
    %create preliminary normalized transients to select baseline
    clear S
    for j=1:size(TRANSD,2)
        S(:,j)=TRANSD(:,j)-mean(TRANSD(:,j));
        S(:,j)=S(:,j)/max(S(:,j));
    end
    plot(S);hold on
    if exist('ST','var')==1
        plot(ST,'k');
    end
    %define baseline for voltage data
    fprintf('SELECT BASELINE, ESC TO ACCEPT.\n')
    X=[];
    button=1;
    while button~=27
        [xi,yi,button]=ginput(1);
            if button~=27
                X=[X;round(xi)];
            end
    end
    %define data sets for fitting and mark in plot
    FITSCANS=[];INTERVALX=[];
    for j=1:2:length(X)-1
        INTERVALX=(X(j):1:X(j+1))';  
        FITSCANS=[FITSCANS;INTERVALX];
        plot(INTERVALX,S(INTERVALX,:),'m')
    end
        
    %define upstroke range
    fprintf('SELECT UPSTROKE RANGE (2 CLICKS).\n')
    [ax,ay]=ginput(2);
    %convert ax into scans
    UPSTROKERANGE=[floor(ax(1)):ceil(ax(2))]';
    %mark upstrokerange
    plot(UPSTROKERANGE,S(UPSTROKERANGE,:),'k')
    %% determine baseline and peak
    UPSTROKE=TRANSD(UPSTROKERANGE,:);
    clear MINSCANS MAXSCANS BASELINESCANS PEAKSCANS BASELINE PEAK
    for j=1:size(UPSTROKE,2)
        T=[UPSTROKERANGE,UPSTROKE(:,j)];
        SORTT=sortrows(T,2);
        MINSCANS(:,j)=SORTT(1:2*escans,1);
        MAXSCANS(:,j)=SORTT(end-escans-1:end,1);
        B=[min(MINSCANS(:,j)):max(MINSCANS(:,j))]';%approximate baseline
        %take only those B scans that are before the maximum
        [val,idx]=max(TRANSD(UPSTROKERANGE,j));idx=idx+UPSTROKERANGE(1)-1;
        BB=[];
        for k=1:length(B)
            if B(k)<idx
                BB=[BB;B(k)];
            end
        end
        if length(BB)>4
            BASELINESCANS(:,j)=[round(mean(BB))-2*escans:round(mean(BB))];
        else
            BASELINESCANS(:,j)=[UPSTROKERANGE(1)-2*escans:UPSTROKERANGE(1)];
        end
        P=[min(MAXSCANS(:,j)):max(MAXSCANS(:,j))]';%approximate maximum
        PEAKSCANS(:,j)=[floor(round(mean(P))-escans/2):ceil(round(mean(P))+escans/2)]';
        BASELINE(j)=mean(TRANSD(BASELINESCANS(:,j),j));
        PEAK(j)=mean(TRANSD(PEAKSCANS(:,j),j));
    end
    %% fit line for voltage measurement
    %perform linear fit
    clear NTRANSD
    for j=1:size(TRANSD,2)
        P=polyfit(FITSCANS,TRANSD(FITSCANS,j),1);
        %draw fit line 
        T=polyval(P,1:size(TRANSD,1))';
        %correct transient
        NTRANSD(:,j)=TRANSD(:,j)-T;
        %correct baseline of transient using defined baseline
        NTRANSD(:,j)=NTRANSD(:,j)-mean(NTRANSD(BASELINESCANS(:,j),j));
        %normalize transient using defined peak scans
        NTRANSD(:,j)=NTRANSD(:,j)/mean(NTRANSD(PEAKSCANS(:,j),j));
    end
    %plot all normalized transients in one graph
    figure(fig);hold off
    TRANST=[1:size(TRANSD,1)]'/scanrate*1.0e3; 
    plot(TRANST,NTRANSD,'k');hold on
    plot(TRANST,0,'r');
    plot(UPSTROKERANGE/scanrate*1e3,NTRANSD(UPSTROKERANGE,:),'m')
    %plot(TRANST,threshold,'r');
    %% measure APD and prepare transients
    APD=zeros(1,size(NTRANSD,2));
    VMAX=zeros(1,size(NTRANSD,2));

    for j=1:size(TRANSD,2)
        [upidx,downidx,T]=apd(NTRANSD(:,j),threshold,UPSTROKERANGE,BASELINESCANS(:,j),scanrate);
        [D0,D1,D2]=smoothdiff(T(UPSTROKERANGE),golaywindow);
        if upidx>0 && downidx>0
            APD(j)=(downidx-upidx)/scanrate*1e3;
        else
            APD(j)=0;
        end
        VMAX(j)=max(D1)*scanrate;
        ALLAPD(i,j)=APD(j);
        ALLVMAX(i,j)=VMAX(j);
        ADDT=[zeros(addscans,1);T;zeros(addscans,1)];
        TI(:,j)=ADDT(floor(upidx+addscans-1)-bscans:floor(upidx+addscans-1)-bscans+totalscans-1);
        if exist('ASTIM','var')==1
            SI(:,j)=ASTIM(floor(upidx+addscans-1)-bscans:floor(upidx+addscans-1)-bscans+totalscans);
        end
        fprintf(['File:',num2str(RAWDATA(i).FILE),' ,#',num2str(j),' APD:',num2str(APD(j)),' ms VMAX:',num2str(VMAX(j)),' 1/s\n']);
        plot(upidx/scanrate*1.0e3,threshold,'r+');
        plot(downidx/scanrate*1.0e3,threshold,'r+');
    end
    
    APDATA(i).RDATA=RAWDATA(i).FDATA;
    APDATA(i).NDATA=NTRANSD;
    APDATA(i).TRANS=TI;
    if exist('SI')==1
        APDATA(i).TSTIM=SI;
        APDATA(i).RSTIM=RAWDATA(i).STIMULUS;
        APDATA(i).ratebpm=RATE(i);
    end
    APDATA(i).threshold=threshold;
    APDATA(i).file=RAWDATA(i).FILE;
    APDATA(i).path=RAWDATA(i).PATH;
    APDATA(i).scanrate=scanrate;
    APDATA(i).APDms=APD;
    askrepeat=input('[RETURN]=accept, [1]=ignore, [2]=repeat >');
    if isempty(askrepeat)==1
        repeat=0;
    else
        if askrepeat==1
            DELDATA=[DELDATA;i];
            repeat=0;
        end
    end 
    end
    %delte entries in APDATA
    APDATA(DELDATA)=[];
end
%% save results
if exist('textpath')==0
    path='C:\';
else
    path=textpath;
end
datapath=uigetdir(path,'DATA FOLDER'); 
datapath=[datapath,'\'];
DESC=input(['Name of data:>'],'s');
apfile=[DESC,'-APDATA','.mat'];
apfiletxt=[DESC,'-APDs','.txt'];
vmaxfiletxt=[DESC,'-VMAXs','.txt'];
apfiletxtrate=[DESC,'-ratesBPM','.txt'];

% LINE=file, COLUMN=regioin
save([datapath,apfile],'APDATA');
save([datapath,apfiletxt],'ALLAPD','-ascii');
save([datapath,vmaxfiletxt],'ALLVMAX','-ascii');
save([datapath,apfiletxtrate],'RATE','-ascii');