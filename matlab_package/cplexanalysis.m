%process cardioplex data
%% User parameter
camrate=input('Camera frame rate (Hz):');
ratiorate=input('Number of frames per ratio:');
importscanrate=camrate/ratiorate;
scanrate=250;%[Hz] >=importscanrate to interpolate data
translength=1250;%[ms] transient length
beforemax=450;% [ms] transient length before upstroke
extremelength=20;% [ms] maximum and minimum length
threshold=0.2;%amplitude threshold of normalized transient
tscans=round(translength*1.0e-3*scanrate)+1;
totalscans=tscans*round(scanrate/importscanrate);%number of scans for the interpolated peak
escans=round(extremelength*1.0e-3*scanrate)+1;
bscans=round(beforemax*1.0e-3*scanrate)+1;
getfrequencies=0;%1=calculate stimulation frequencies
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

if 0<stimcol && stimcol<=2 && getfrequencies==1
    %go through data and extract stimulus and frequency from each entry
    hdl = waitbar(0,['CALCULATING STIMULATION FREQUENCIES']);
    R=zeros(length(RAWDATA),1);
    for i=1:length(RAWDATA)
        if stimcol==1
            RAWDATA(i).STIMULUS=RAWDATA(i).DATA(:,1);
            RAWDATA(i).DATA(:,1)=[];
        else
            RAWDATA(i).STIMULUS=RAWDATA(i).DATA(:,end);
            RAWDATA(i).DATA(:,end)=[];
        end
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
    Hopts.Fs=importscanrate;
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
    R(i)=maxf;
    waitbar(i/length(RAWDATA));
    end
    close(hdl);
end 
fprintf('Stimulation frequencies:\n');
R
%% select peaks to analyze
fig=figure('NumberTitle','off');
PEAKDATA=[];DELDATASET=[];
fprintf('SELECT UPSTROKERANGE (2 CLICKS EACH TRANSIENT).\n');
for i=1:length(RAWDATA)
    repeat=1;
    while repeat==1
    ALLCOL=RAWDATA(i).DATA;%data in all columns of set i
    ALLCOLPEAKS=zeros(totalscans,size(ALLCOL,2));
    ALLCOLTIME=zeros(size(ALLCOLPEAKS));
    ALLDIAST=zeros(size(ALLCOL,2),1);
    ALLSYST=zeros(size(ALLCOL,2),1);
    ALLAMP=zeros(size(ALLCOL,2),1);
    ALLDUR=zeros(size(ALLCOL,2),1);
    ALEFTTIME=zeros(size(ALLCOL,2),1);
    ARIGHTTIME=zeros(size(ALLCOL,2),1);
    
    if stimcol>0
        ALLSTIM=zeros(size(ALLCOLPEAKS));
    end

    for j=1:size(ALLCOL,2)
      TRANSD=ALLCOL(:,j);%data
      TRANST=[1:length(TRANSD)]'/importscanrate*1.0e3;
      
      %only for field stimulation can one assume simultaneous excitation
      %(remove the if condition for propagated transients)
      if j==1 
        figure(fig);
        set(fig,'Name',['UPSTROKE RANGE FILE: ',num2str(RAWDATA(i).FILE),' COL ',num2str(j),' OF ',num2str(size(ALLCOL,2))]);
        plot(TRANST,TRANSD,'k'); hold on
        [ax,ay]=ginput(2);
        %convert ax into scans
        X=round([ax(1),ax(2)]*1.0e-3*importscanrate);
        UPSTROKERANGE=[X(1)-1:X(2)+1]';
        %mark upstrokerange
        plot(TRANST(UPSTROKERANGE),TRANSD(UPSTROKERANGE),'b')
      end
      
      %measurements:sort UPSTROKERANGE
      UPSTROKE=[UPSTROKERANGE,TRANSD(UPSTROKERANGE)];
      SORTUPSTROKE=sortrows(UPSTROKE,2);

      %get range of scans with smallest amplitudes
      MINSCANS=SORTUPSTROKE(1:escans,:);%get the scans with the smallest amplitudes
      MINUPSTROKERANGE=[min(MINSCANS(:,1)):max(MINSCANS(:,1))];
      MINUPSTROKE=[MINUPSTROKERANGE',TRANSD(MINUPSTROKERANGE)];
      
      %get range of scans with largest amplitudes
      MAXSCANS=SORTUPSTROKE(end-escans+1:end,:);%get scans with largest amplitudes
      MAXUPSTROKERANGE=[min(MAXSCANS(:,1)):max(MAXSCANS(:,1))];
      MAXUPSTROKE=[MAXUPSTROKERANGE',TRANSD(MAXUPSTROKERANGE)];
           
      diastlevel=mean(MINUPSTROKE(:,2));
      systlevel=mean(MAXUPSTROKE(:,2));
      amplevel=systlevel-diastlevel;
      
      %normalize transient and measure duration
      N=(TRANSD-diastlevel)/amplevel;
      
      %left index
        idx=MINUPSTROKE(1,1);
        while N(idx)<threshold;idx=idx+1;end;left=idx-1;
        %refine left index
        XDATA=[idx-2:idx+2]';
        FITXDATA=[XDATA(1):0.1:XDATA(end)]';
        P=polyfit(XDATA,N(XDATA),1);
        FITYDATA=polyval(P,FITXDATA);
        fitidx=1;while FITYDATA(fitidx)<threshold;fitidx=fitidx+1;end;fineleft=FITXDATA(fitidx-1);
      
      %right index
        idx=MAXUPSTROKE(end,1);
        while N(idx)>threshold && idx<length(N);idx=idx+1;end;right=idx-1;
        if idx==length(N)
            fineright=fineleft;
        else
            XDATA=[idx-20:idx+1]';
            FITXDATA=[XDATA(1):0.1:XDATA(end)]';
            P=polyfit(XDATA,N(XDATA),1);
            FITYDATA=polyval(P,FITXDATA);
            fitidx=2;
            while fitidx<length(FITYDATA)-1 ...
                && FITYDATA(fitidx)<FITYDATA(fitidx-1) ...
                && FITYDATA(fitidx)>threshold;
                fitidx=fitidx+1;
            end;
            if 2<fitidx && fitidx<length(FITYDATA)
                fineright=FITXDATA(fitidx+1);
            else
                fineright=fineleft;
            end
        end
        
                    
      %mark measurements
      lefttime=fineleft/importscanrate*1.0e3;
      righttime=fineright/importscanrate*1.0e3;
      dlevel=diastlevel+amplevel*threshold;
      if j==1
        figure(fig)
        plot(TRANST(MINUPSTROKE(:,1)),MINUPSTROKE(:,2),'g+');
        plot(TRANST,diastlevel,'g');
        plot(TRANST(MAXUPSTROKE(:,1)),MAXUPSTROKE(:,2),'r+');
        plot(TRANST,systlevel,'r');
        plot(lefttime,dlevel,'r+')
        plot(righttime,dlevel,'r+')
        plot([lefttime:righttime]',dlevel,'r');hold off;
      end
      %cut peak and stimulus
      startleft=round(mean(MAXUPSTROKE(:,1)-bscans))-1;
      if startleft<1
          TRANSDADD=[mean(TRANSD(1:5))*ones(abs(startleft+1),1);TRANSD];
          CUTSCANS=[1:tscans+1]';
          CUTPEAK=TRANSDADD(CUTSCANS);
      else
          CUTSCANS=[startleft:startleft+tscans]';
          CUTPEAK=TRANSD(CUTSCANS);
      end
      TIME=([1:length(CUTPEAK)]'-bscans)*1.0e3/importscanrate;
      ITIME=[TIME(1):1/scanrate*1.0e3:TIME(end)]';
      ITIME=ITIME(1:totalscans);
      ICUTPEAK=interp1(TIME,CUTPEAK,ITIME,'cubic');
      ICUTPEAK=ICUTPEAK(1:totalscans,1);
      %cut stimulus
      if stimcol>0
          STIM=RAWDATA(i).STIMULUS;
          CUTSTIM=STIM(CUTSCANS);
          ICUTSTIM=interp1(TIME,CUTSTIM,ITIME,'cubic');
          ICUTSTIM=ICUTSTIM(1:totalscans,1);
      end
      %save data
      ALLCOLPEAKS(:,j)=ICUTPEAK;
      ALLSTIM(:,j)=ICUTSTIM;
      ALLCOLTIME(:,j)=ITIME;
      ALLDIAST(j)=diastlevel;
      ALLSYST(j)=systlevel;
      ALLAMP(j)=amplevel;
      ALLDUR(j)=righttime-lefttime;
      ALEFTTIME(j)=lefttime;
      ARIGHTTIME(j)=righttime;
      
      if j==1
        set(fig,'Name','RETURN TO CONTINUE');
        fprintf('RETURN TO CONTINUE\n');
        pause
      end
    end
    PEAKDATA(i).DATA=ALLCOLPEAKS;
    PEAKDATA(i).STIM=ALLSTIM;
    PEAKDATA(i).TIME=ITIME;
    PEAKDATA(i).rate=RAWDATA(i).rate;
    PEAKDATA(i).DIAST=ALLDIAST;
    PEAKDATA(i).SYST=ALLSYST;
    PEAKDATA(i).AMP=ALLAMP; 
    PEAKDATA(i).DURMS=ALLDUR;
    PEAKDATA(i).FILE=RAWDATA(i).FILE;
    PEAKDATA(i).PATH=RAWDATA(i).PATH;
    %plot results for all hearts on tile plot (3 x 3)
    tilefig=figure('Units','normalized','Position',[0.1,0.1,0.75,0.75],'Name','Summary');
    set(tilefig,'Name',['Summary file: ',num2str(PEAKDATA(i).FILE)])
    for j=1:size(ALLCOLPEAKS,2)
        subplot(3,3,j)
        plot(PEAKDATA(i).TIME,PEAKDATA(i).DATA(:,j),'k');
        hold on
        plot([PEAKDATA(i).TIME(1):PEAKDATA(i).TIME(end)]',PEAKDATA(i).DIAST(j),'g')
        plot([PEAKDATA(i).TIME(1):PEAKDATA(i).TIME(end)]',PEAKDATA(i).SYST(j),'r')
        hold off
    end
    askaccept=input('[RETURN]=ACCEPT MEASUREMENT (1)=IGNORE MEASUREMENT (2)=TRY AGAIN>');
    if isempty(askaccept)==1
        repeat=0;
    elseif askaccept==1
        DELDATASET=[DELDATASET;i];
        repeat=0;
    end
    end    
end
close(fig);
%% Delete datasets
if isempty(DELDATASET)~=1
    PEAKDATA(DELDATASET)=[];
end
%% Prepare output
%Determine number of data columns for each file
for i=1:size(PEAKDATA,2)
    NUMCOL(i)=size(PEAKDATA(i).DATA,2);
end
%Report missing data
fprintf(['Maximum number of transients per file in data: ',num2str(max(NUMCOL)),'.\nReporting missing data below:\n']);
for i=1:size(PEAKDATA,2)
    if size(PEAKDATA(i).DATA,2)<max(NUMCOL)
        fprintf([num2str(size(PEAKDATA(i).DATA,2)),' transients in file ',num2str(PEAKDATA(i).FILE),'\n\n']);
    end
end

%Prepare Output Matrices
TRC=zeros(length(PEAKDATA(1).DATA),max(NUMCOL));
T=ITIME;
DIA=zeros(length(PEAKDATA),max(NUMCOL));
SYS=zeros(size(DIA));
AMP=zeros(size(DIA));
DURMS=zeros(size(DIA));
RATE=zeros(length(PEAKDATA),1);

%Fill Output Matrices
%REGIONS ON EACH HEART

for j=1:max(NUMCOL)
%HEARTS
TEMP=[];
    for i=1:length(PEAKDATA)
        RATE(i)=PEAKDATA(i).rate;
        if j<=size(PEAKDATA(i).DATA,2)
            TEMP=[TEMP,PEAKDATA(i).DATA(:,j)];
            DIA(i,j)=PEAKDATA(i).DIAST(j);
            SYS(i,j)=PEAKDATA(i).SYST(j);
            AMP(i,j)=PEAKDATA(i).AMP(j);
            DURMS(i,j)=PEAKDATA(i).DURMS(j);
        end
    end
    %AVERAGE TRANSIENTS IN REGION j
    TRC(:,j)=mean(TEMP,2);
end
%% save data
TRCFILE='MEAN-TRANSIENTS.txt';
TIMEFILE='TIME.txt';
DIAFILE='DIASTOLIC-LEVELS.txt';
SYSFILE='SYSTOLIC-LEVELS.txt';
AMPFILE='RELEASE-LEVELS.txt';
DURFILE=['TRANSIENTDURATION',num2str(threshold*100),'-MS.txt'];
RATEFILE='STIMFREQUENCIES.txt';

save([PEAKDATA(1).PATH,TRCFILE],'TRC','-ascii','-tabs');
save([PEAKDATA(1).PATH,TIMEFILE],'T','-ascii','-tabs');
save([PEAKDATA(1).PATH,DIAFILE],'DIA','-ascii','-tabs');
save([PEAKDATA(1).PATH,SYSFILE],'SYS','-ascii','-tabs');
save([PEAKDATA(1).PATH,AMPFILE],'AMP','-ascii','-tabs');
save([PEAKDATA(1).PATH,DURFILE],'DURMS','-ascii','-tabs');
save([PEAKDATA(1).PATH,RATEFILE],'RATE','-ascii','-tabs');

%print file names
fprintf('FILE NAMES:\n')
for i=1:length(PEAKDATA)
    fprintf(['HEART (ROW):',num2str(i),' NAME:',num2str(PEAKDATA(i).FILE),' \n']);
end
    

fprintf('Description of output files:\n');
fprintf('MEAN-TRANSIENTS.txt: Mean transient across all files sorted by columns. \n');
fprintf('TIME.txt: Time for MEAN-TRANSIENTS.txt in ms. Maximum at 0.\n');
fprintf('DIASTOLIC-LEVELS.txt: Lowest-ratios matrix (line,column)=(file,transient). Zero for missing transients.\n');
fprintf('SYSTOLIC-LEVELS.txt: Highest-ratios matrix (line,column)=(file,transient). Zero for missing transients.\n');
fprintf('RELEASE-LEVELS.txt: (Highest-Lowest)-ratios matrix (line,column)=(file,transient). Zero for missing transients.\n');
fprintf('TRANSIENTDURATIONXX-MS.txt: Transient duration at XX level in milliseconds\n');