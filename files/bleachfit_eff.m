%clear RAWDATA NRAWDATA
%% pick sample pixels
[FBRIDX,map]=gray2ind(FR,256);
RGBFBR=ind2rgb(FBRIDX,gray(256));
if exist('PIXELCOORDINATES')==0
    [PIXELCOORDINATES,RGBFBR]=pixelselect(FR,[],[],FILTERDATA);
end
%% select sample pixels
FITPIXELS=zeros(size(FILTERDATA,3),size(PIXELCOORDINATES,1));
for i=1:size(PIXELCOORDINATES,1)
    FITPIXELS(:,i)=-squeeze(FILTERDATA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
end
%correct offset, normalize and display sample pixels
NFITPIXELS=zeros(size(FITPIXELS));%normalized sample pixels
for i=1:size(FITPIXELS,2)
        NFITPIXELS(:,i)=FITPIXELS(:,i)-mean(FITPIXELS(:,i));
        NFITPIXELS(:,i)=NFITPIXELS(:,i)/max(NFITPIXELS(:,i));
end
%plot sample time intervals
alldata=figure('Name', 'SELECT EPISODE, 2 CLICKS');
plot(NFITPIXELS,'k');hold on
if exist('STIM')==1 && isempty(STIM)==0
    %plot(STIM,'b')
end
%plot previous interval selection
if exist('EPISODEX')==1 && isempty(EPISODEX)==0
    plot(EPISODEX,NFITPIXELS(EPISODEX,:),'r');
end
%[ax,ay]=ginput(2);

%hard code to select entire length for analysis for efficiency:
ax(1)=ceil(.001*length(NFITPIXELS));
ax(2)=ceil(.999*length(NFITPIXELS));

EPISODEINTERVAL=round([ax(1),ax(2)]);
EPISODEFRAME=[EPISODEINTERVAL(1),EPISODEINTERVAL(2)]+firstframenumber-1;
EPISODEX=[EPISODEINTERVAL(1):EPISODEINTERVAL(2)];
figure(alldata);
plot(EPISODEX,NFITPIXELS(EPISODEX,:),'r');
EPISODE=FILTERDATA(:,:,EPISODEX);
EPISODETRACE=NFITPIXELS(EPISODEX,:);
if exist('STIM')==1 && isempty(STIM)==0
    EPISODESTIM=STIM(EPISODEX,1);
end
%% choose intervals for fitting

%OBTAIN FITTING INTERVALS
hdimage=figure('Name','FITTING INTERVALS');
plot(EPISODETRACE(:,1),'b');hold on
plot(EPISODETRACE(:,2),'r')

%fitmethod=input('1=select baseline 2=select upstrokeranges>');
%hard code to "1" for efficiency:

fitmethod=1

if fitmethod==1    
    button=1;
    if exist('XPOINTS','var')==1 && isempty(XPOINTS)==0
        %ask to enter new baseline if X exists
        enterbaseline=input('keep baseline [return]:');
        if isempty(enterbaseline)==1
            enterbaseline=0;
        else
            enterbaseline=1;
        end
    else
        enterbaseline=1;
    end
    
    %enter new baseline
    if enterbaseline==1    
        XPOINTS=[];
        while button~=27
            [xi,yi,button]=ginput(1);
            if button~=27
                XPOINTS=[XPOINTS;round(xi)];
            end
        end
    end
    
    
    saveas(gcf,firstfile(1:length(firstfile)-4),'tiff')
        
    %define data sets for fitting and mark in plot
    XDATA=[];INTERVALX=[];
    figure(hdimage);hold on
    for i=1:2:length(XPOINTS)-1
        INTERVALX=(XPOINTS(i):1:XPOINTS(i+1))';  
        XDATA=[XDATA;INTERVALX];
        plot(INTERVALX,EPISODETRACE(INTERVALX,:),'m')
    end

else
    UPFIT=[];
    for i=1:2
        fprintf(['select range ',num2str(i),': baseline + upstrokes\n']);
        [ax,ay]=ginput(2);
        UPFIT(i).R=[round(ax(1)):round(ax(2))]';
        X=UPFIT(i).R;Y=EPISODETRACE(X,:);
        plot(X,Y,'k');
    end
end

%% fit all pixels
BLEACHDATA=zeros(size(EPISODE));
OPTIONS=optimset('display','off','MaxFunEvals',200,'MaxIter',200,'TolFun',1e-6,'TolX',1e-6);
maxvalue=65536;%maximum data value is 2^16
P(2)=100*size(BLEACHDATA,3);
btime=75;%time corresponding to baselinescans [ms]
bscans=btime*scanrate*1.0e-3;
if fitmethod==1
    BASELINESCANS=zeros(size(EPISODE,1),size(EPISODE,2),length(XDATA));
else
    BASELINESCANS=zeros(size(EPISODE,1),size(EPISODE,2),bscans*length(UPFIT));%Baseline scans used in fit
end
apixel=ceil(10*scanrate/2000);
hdl = waitbar(0,['FITTING BASELINE FOR ALL PIXELS']);
for j=1:size(EPISODE,2)
    parfor i=1:size(EPISODE,1)
        PIXEL=squeeze(EPISODE(i,j,:));
%         if fitmethod==2
%             %find XDATA for each pixel
%             UPSTROKEFIT=struct(UPFIT);
%             XDATA=[];
%             for k=1:length(UPSTROKEFIT)
%                 R=UPSTROKEFIT(k).R;%get upstroke range k
%                 UPSTROKEAMP=[];
%                 for m=R(1):R(end);
%                     UPSTROKEAMP=[UPSTROKEAMP;[m,(PIXEL(m))]];
%                 end
%                 %sort amplitudes
%                 SORTUPSTROKEAMP=sortrows(UPSTROKEAMP,2);
%                 %get bscans
%                 XDATA=[XDATA;SORTUPSTROKEAMP(end-bscans+1:end,1)];
%             end
%             XDATA=sort(XDATA);
%         end
    BASELINESCANS(i,j,:)=XDATA;
    %fit PIXEL
%% Baseline fit
METHOD='linear';   
FITBASELINE=baselinefit(PIXEL,XDATA,METHOD);
%method is "linear","poly" or "exp"
NEWPIXEL=(PIXEL-FITBASELINE);
BLEACHDATA(i,j,:)=NEWPIXEL;%replace pixel in new data matrix
    end
    waitbar(j/size(EPISODE,2))
end

%show fit for example pixels

close(hdl);
hdimage=figure('Name','CLICK 4 PIXELS TO VIEW RESULT');
image(RGBFBR)

%choose 4 pixels for plotting
VPIXEL=[];Y=[];
for i=1:4
 [xi,yi,enter]=ginput(1);
 VPIXEL(i,:)=[round(yi),round(xi)];
end
figure;hold on
for i=1:size(VPIXEL,1)
subplot(2,2,i);hold on
TRACE=squeeze(BLEACHDATA(VPIXEL(i,1),VPIXEL(i,2),:));
plot(TRACE);
%mark baseline in figure
%plot(1:length(TRACE),mean(TRACE(squeeze(BASELINESCANS(VPIXEL(i,1),VPIXEL(i,2),1:bscans)))),'r')
hold off
end

%prepare corrected EPISODETRACE
for i=1:size(PIXELCOORDINATES,1)
    XDATA=squeeze(BASELINESCANS(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
    T=squeeze(BLEACHDATA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
    T=T-mean(T(XDATA));
    EPISODETRACE(:,i)=T/min(T);
end
%% save data to maximize memory for NORMDATA array
if exist('datapath')==0 
    datapath=stackpath;
else
    if length(datapath)<2
        datapath=stackpath;
    end
end

%datapath=uigetdir(datapath,'SELECT FOLDER FOR IMAGES');
%datapath=[datapath,'\'];
% imagepath will be set equal to datapath so the images are saved to the same folder for efficiency


[fpath, fname, fext] = fileparts(stackfile);
bleachfile=[stackfile(1:end-length(fext)),'-B.mat'];
bleachpathfile=[datapath,bleachfile];


fprintf(['SAVING DATA \n'])
bleachpathfile
if exist('EPISODESTIM')==0
    EPISODESTIM=[]; 
end
SAVESTIM=STIM;%copy original stimulus 
STIM=EPISODESTIM;%replace STIM for saving

save(bleachpathfile,'FILTERDATA','STIM','EPISODESTIM','BLEACHDATA','BASELINESCANS',...
    'EPISODEFRAME','EPISODEX','EPISODETRACE','scanrate','stackfile','datapath',...
    'bleachfile','firstframenumber','lastframenumber','FITPIXELS','PIXELCOORDINATES',...
    'FR','RGBFBR','stackpath','stackfile','fitmethod','lpcutoff','spatpix','medianfilter');
STIM=SAVESTIM;clear SAVESTIM%restore original stimulus