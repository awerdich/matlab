%divide data in episodes and fit baseline
%% pick sample pixels
if exist('PIXELCOORDINATES')==0
    %[PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,SCATTER,ASIGNALPIXELS,DATA)
    [PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,[],[],FR);
end
%% select sample pixels
FITPIXELS=zeros(size(FR,3),size(PIXELCOORDINATES,1));
for i=1:size(PIXELCOORDINATES,1)
    FITPIXELS(:,i)=squeeze(FR(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
end
%correct offset, normalize and display sample pixels
NFITPIXELS=zeros(size(FITPIXELS));%normalized sample pixels
for i=1:size(FITPIXELS,2)
        NFITPIXELS(:,i)=FITPIXELS(:,i)-mean(FITPIXELS(:,i));
        NFITPIXELS(:,i)=NFITPIXELS(:,i)/max(NFITPIXELS(:,i));
end
% add stimulus
NFITPIXELS(:,end+1) = STIM;
%% select episode
%plot selected traces
alldata=figure('Name', 'SELECT EPISODE');
plot(NFITPIXELS);hold on
%plot previous interval selection
if exist('EPISODEX')==1 && isempty(EPISODEX)==0
    plot(EPISODEX,NFITPIXELS(EPISODEX,:),'r');
end
[ax,ay]=ginput(2);
%ax(1)=1;ax(2)=size(NFITPIXELS,1);
EPISODEINTERVAL=round([ax(1),ax(2)]);
EPISODEFRAME=[EPISODEINTERVAL(1),EPISODEINTERVAL(2)];
EPISODEX=[EPISODEINTERVAL(1):EPISODEINTERVAL(2)]';
figure(alldata);
plot(EPISODEX,NFITPIXELS(EPISODEX,:),'r');
EPISODETRACE=NFITPIXELS(EPISODEX,:);
if exist('STIM')==1 && isempty(STIM)==0
    EPISODESTIM=STIM(EPISODEX,1);
end
%% define EPISODE based on frame selection
EPISODEW1=FW1(:,:,EPISODEX);
EPISODEW2=FW2(:,:,EPISODEX);
EPISODER=FR(:,:,EPISODEX);
if exist('FCA')==1
    EPISODECA=FCA(:,:,EPISODEX);
end
%% Correct baseline drift of episode
%askfit=input('Correct baseline drift? (1=YES)>');
askfit=0;
if isempty(askfit)==0 && askfit==1
    figure;
    plot(EPISODETRACE(:,1),'b');hold on
    plot(EPISODETRACE(:,2),'r')

    UPSTROKERANGE=[];
    for i=1:2
    fprintf(['select range ',num2str(i),': baseline + upstrokes\n']);
    [ax,ay]=ginput(2);
    UPSTROKERANGE(i).R=[round(ax(1)):round(ax(2))]';
    X=UPSTROKERANGE(i).R;Y=EPISODETRACE(X,:);
    plot(X,Y,'k');
    end

    %FIT ROUTINE
    EPISODEFITR=zeros(size(EPISODER));
    if exist('FCA')==1
    EPISODEFITCA=zeros(size(EPISODECA));
    end
    baselinescans=10;
    BASELINESCANS=zeros(size(EPISODER,1),size(EPISODER,2),baselinescans*length(UPSTROKERANGE));%BASELINESCANS FOR ALL PIXELS
    apixel=2;
    hdl = waitbar(0,['FITTING DATA']);
    for j=1:size(EPISODER,2)
        for i=1:size(EPISODER,1)
    %% FIT PIXEL
            PIXEL=squeeze(EPISODER(i,j,:));
            if exist('FCA')==1
                PIXEL=[PIXEL;squeeze(EPISODECA(i,j,:))];
            end
            T=PIXEL(:,1);
            %define baselines of T
            BASE=zeros(baselinescans,length(UPSTROKERANGE));%BASELINESCANS FOR THIS PIXEL
            XDATA=zeros(baselinescans*length(UPSTROKERANGE),1);
            for k=1:length(UPSTROKERANGE)
                R=UPSTROKERANGE(k).R;%get upstroke range k
                UPSTROKEAMP=[];
                for m=R(1):R(end);
                    UPSTROKEAMP=[UPSTROKEAMP;[m,mean(T(m-round(apixel/2):m+round(apixel/2)-1))]];
                end
                %sort amplitudes
                SORTUPSTROKEAMP=sortrows(UPSTROKEAMP,2);
                %get baselinescans
                BASE(:,k)=SORTUPSTROKEAMP(1:baselinescans,1);
                XDATA(baselinescans*(k-1)+1:baselinescans*k,1)=BASE(:,k);
            end
            XDATA=sort(XDATA);
            BASELINESCANS(i,j,:)=XDATA;
            %fit PIXEL
            FITPIXEL=zeros(size(PIXEL));%fitted PIXEL
            for k=1:size(PIXEL,2)
                P(k,:)=polyfit(XDATA,PIXEL(XDATA,k),1);
                %USE FIRST BASELINE
                %FITPIXEL(:,k)=PIXEL(:,k)-((polyval(P(k,:),1:size(PIXEL,1))')-mean(PIXEL(XDATA(1:baselinescans),k),1));
                %USE SECOND BASELINE
                FITPIXEL(:,k)=PIXEL(:,k)-((polyval(P(k,:),1:size(PIXEL,1))')-mean(PIXEL(XDATA(baselinescans+1:end),k),1));
            end
            %write FITPIXEL in EPISODEFIT
            EPISODEFITR(i,j,:)=FITPIXEL(:,1);
            if exist('FCA')==1
            EPISODEFITCA(i,j,:)=FITPIXEL(:,2);
            end
        end
        waitbar(j/size(EPISODER,2));
    end
    close(hdl);
    %display results
    BASE=zeros(baselinescans*length(UPSTROKERANGE),size(PIXELCOORDINATES,1));
    for i=1:size(PIXELCOORDINATES,1)
        %FITPIXELS(:,i)=squeeze(FCA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
        PIXEL(:,i)=squeeze(EPISODER(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));%PRE-FIT PIXELS
        FITPIXEL(:,i)=squeeze(EPISODEFITR(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));%POST-FIT PIXELS
        BASE(:,i)=squeeze(BASELINESCANS(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
    end
    %plot results for PIXELCOORDINATES
    figure
    plot(PIXEL(:,1),'b');hold on
    plot(FITPIXEL(:,1),'b');
    plot(BASE(:,1),FITPIXEL(BASE(:,1),1),'b+');
    plot(1:size(FITPIXEL,1),mean(FITPIXEL(BASE(1:baselinescans,1),1)),'b');

    plot(PIXEL(:,2),'r');
    plot(FITPIXEL(:,2),'r');
    plot(BASE(:,2),FITPIXEL(BASE(:,2),2),'r+');
    plot(1:size(FITPIXEL,1),mean(FITPIXEL(BASE(1:baselinescans,2),2)),'r');
    
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
    TRACE1=squeeze(EPISODER(VPIXEL(i,1),VPIXEL(i,2),:));
    BASE1=squeeze(BASELINESCANS(VPIXEL(i,1),VPIXEL(i,2),:));
    TRACE2=squeeze(EPISODEFITR(VPIXEL(i,1),VPIXEL(i,2),:));
    BASE2=squeeze(BASELINESCANS(VPIXEL(i,1),VPIXEL(i,2),:));

    plot(TRACE1,'k');
    plot(BASE1,TRACE1(BASE1),'k+');
    plot(TRACE2,'r');
    plot(BASE2,TRACE2(BASE2),'r+');
    plot(1:length(TRACE2),mean(TRACE2(BASE2)),'m');

    hold off
    end 
else
    askfit=0;
end
%% save data to maximize memory for NORMDATA array

if exist('datapath')==0 || isempty(datapath)==1 || length(datapath)<4
    path=stackpath;
    datapath=uigetdir(path,'SELECT FOLDER TO SAVE DATA FILES');
    datapath=[datapath,'\'];
end
[fpath, fname, fext] = fileparts(stackfile);
bleachfile=[stackfile(1:end-length(fext)),'-B.mat'];
bleachpathfile=[datapath,bleachfile];

%save data after fitting
DATAST(1).spatpix=spatpix;
DATAST(1).lpcutoff=lpcutoff;
DATAST(1).FW1=FW1;%
DATAST(1).FW2=FW2;%
DATAST(1).STIM=STIM;
DATAST(1).EPISODESTIM=EPISODESTIM;
DATAST(1).EPISODEW1=EPISODEW1;
DATAST(1).EPISODEW2=EPISODEW2;
DATAST(1).EPISODER=EPISODER;

save(bleachpathfile,'DATAST');