%% set up new database file
if exist('FDATABASE','var')==0
    fprintf('NEW DATABASE GENERATED \n')
    h = msgbox('No Database in Memory','WARNING')
    FDATABASE=[];id=1;
else
    id=size(FDATABASE,2);
    size(FDATABASE)
    fprintf(['previous heart, id:',num2str(id),' name:',num2str(FDATABASE(id).name),'\n']);
    nextheart=input('This heart number:');
    if isempty(nextheart)==0
        id=nextheart;
    else
        return
    end
end
%% general information
if exist('magnification','var')==0;mapconfig;end
%ask for data entry
editinfo=0;
%cases in which we need to edit heart info
if isfield(FDATABASE,'age')==0 
    editinfo=1;
else
    if length(FDATABASE)<id || isempty(FDATABASE(id).age)==1 || isempty(FDATABASE(id).date)==1
        editinfo=1;
    else
        FDATABASE(id)
        askedit=input('Edit heart info <1>=YES <RETURN>=NO :');
        if askedit==1
            editinfo=1;
        end
    end
end
if editinfo==1
    FDATABASE(id).age=input('age (hpf):');
    FDATABASE(id).date=input('date (yymmdd):');
    FDATABASE(id).heart=input('heart number:');
    FDATABASE(id).rate=input('HR (bpm):');
    FDATABASE(id).sinusloc=input('impulse initiation 0=atrium 1=av 2=ventricle 3=fstim:');
    if FDATABASE(id).sinusloc==3
        FDATABASE(id).stim=FDATABASE(id).rate;
    else
        FDATABASE(id).stim=input('stimulation rate 0:no stimulation (bpm):');
    end
    
    %stretch information
%     FDATABASE(id).stretchregion=input('Area stretched (1) Atrium, (2) Ventricle:');
%     FDATABASE(id).stretchtime=input('stretch time (minutes):');
%     FDATABASE(id).recovertime=input('recovery time (minutes):');

%    FDATABASE(id).comment=input('Drug(Name)-Concentration(uM):','s');
%    FDATABASE(id).path=datapath;
%    FDATABASE(id).camrate=camrate;
%    FDATABASE(id).fratio=fratio;
%    FDATABASE(id).scanrate=scanrate;
%    FDATABASE(id).name=stackfile;
%    FDATABASE(id).pixelcalfactorx=pixelcalfactor_x;
%    FDATABASE(id).pixelcalfactory=pixelcalfactor_y;
%    FDATABASE(id).magnification=magnification;
%    FDATABASE(id).sideport=sideport;
%    FDATABASE(id).lpcutoff=lpcutoff;
%    FDATABASE(id).spatpix=spatpix;
%    FDATABASE(id).softwareversion=pwd;%software version
end
%% save tracedata
fprintf(['Filename:',FDATABASE(id).name,' \n']);
region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
clear regionv
if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
|| strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
|| strcmp(region,'vic')==1 || strcmp(region,'voc')==1

    regionframedata=[region,'_DATA'];%ALL TRANSIENTS IN FRAME (i,j,:)
    regionmeantrace=[region,'_TRACE'];%MEAN FRAME TRANSIENT
    regionframex=[region,'_framex'];%frame width in pixels
    regionframey=[region,'_framey'];%frame height
    regionframeumwidth=[region,'_frameumwidth'];%frameumwidth
    regionframeumheight=[region,'_frameumheight'];%frameumheight
    regionframecenter=[region,'_framecenter'];%frame center locatin
    regionarea=[region,'_area'];%estimate of active tissue area in frame (SIGNALS>0)
    regionlow=[region,'_diast'];%AVERAGE (MINIMUM OF EACH TRANSIENT)
    regionhigh=[region,'_syst'];%AVERAGE (MAXIMUM OF EACH TRANSIENT)
    regionamp=[region,'_amp'];%AVERAGE (AMPLITUDE OF EACH TRANSIENT)
    regionframeij=[region,'_frameij'];%pixel coordinates for each transient used to calculate average
    regiondur=[region,'_durms',num2str(durthreshold*100)];%AVERAGE TRANSIENT DURATION (ms) 
    regionlowmatrix=[region,'_RLOWROI'];%diastolic ROI matrix
    regionhighmatrix=[region,'_RHIGHROI'];%systolic ROI matrix
    regiondurmatrix=[region,'_DUROI'];%duration ROI matrix
else
    fprintf('Wrong input!\n');
    return;
end

FDATABASE(id).(regionframex)=FrameXLim;
FDATABASE(id).(regionframey)=FrameYLim;
FDATABASE(id).(regionframeumwidth)=frameumwidth;
FDATABASE(id).(regionframeumheight)=frameumheight;
FDATABASE(id).(regionframecenter)=FRAMECENTER;

%get mean transient
T=[];A=[];B=[];FRAMEij=[];LOW=[];HIGH=[];AMP=[];U=[];KS=[];DURij=[];
roiheight=FrameYLim(2)-FrameYLim(1)+1;
roiwidth=FrameXLim(2)-FrameXLim(1)+1;
RLOWROI=zeros(roiheight,roiwidth);
RHIGHROI=zeros(size(RLOWROI));
DUROI=zeros(size(RLOWROI));

totalscans=2000;
threshold=0.5;
for i=FrameYLim(1):FrameYLim(2)
    for j=FrameXLim(1):FrameXLim(2)
        %ROI coordinates
        m=i-FrameYLim(1)+1;
        n=j-FrameXLim(1)+1;
        if SIGNAL(i,j)>0
            %save matrix data
            RLOWROI(m,n)=RLOW(i,j);
            RHIGHROI(m,n)=RHIGH(i,j);
            DUROI(m,n)=DUR(i,j);
            %get transient
            if exist('EPISODECA')==1
                A=squeeze(EPISODECA(i,j,:));
            else
                if exist('EPISODEFITR')==1 && askfit==1
                    A=squeeze(EPISODEFITR(i,j,:));
                else
                    A=squeeze(EPISODER(i,j,:));
                end 
            end
            baseline=mean(A(squeeze(BSCANS(i,j,:))));
            %add data to make sure transients are long enough
            A=[ones(round(totalscans/2),1).*mean(A);A;ones(round(totalscans/2),1).*mean(A)];
            if exist('EPISODESTIM','var')==1 && isempty(EPISODESTIM)==0
                ST=[ones(round(totalscans/2),1).*mean(EPISODESTIM);EPISODESTIM;ones(round(totalscans/2),1).*mean(EPISODESTIM)];
            end
            U=round(totalscans/2)+UPSTROKERANGE;
            %determine upstroke after minimum and align
            [mval,midx]=min(A(U));midx=midx+U(1)-1;%minimum
            U2=[midx-10:U(end)]';
            if U2(1)>0
                [val,idx]=max(A(U2));idx=idx+U2(1)-1;%maximum
            else
                [val,idx]=max(A(U));idx=idx+U(1)-1;%maximum
            end
            NA=A-baseline;NA=NA/NA(idx);%normalize transient
            k=idx;while NA(k)>threshold;k=k-1;end;k=k+1;%scan of 92% maximum
            B=A(k-round(totalscans/2)+1:k+round(totalscans/2));
            T=[T,B];
            %save index to determine mean k to cut stimulus            
            KS=[KS;k];
            %get other measurements
            FRAMEij=[FRAMEij;[i,j]];
            if DUR(i,j)>0
                DURij=[DURij;[DUR(i,j),i,j]];
            end
            if exist('CALOW')
                LOW=[LOW;CALOW(i,j)];
                HIGH=[HIGH;CAHIGH(i,j)];
                AMP=[AMP;(CAHIGH(i,j)-CALOW(i,j))];
            else
                LOW=[LOW;RLOW(i,j)];
                HIGH=[HIGH;RHIGH(i,j)];
                AMP=[AMP;(RHIGH(i,j)-RLOW(i,j))];
            end
        end
    end
end

FDATABASE(id).(regionframedata)=T;
FDATABASE(id).(regionmeantrace)=mean(T,2);
FDATABASE(id).(regionarea)=size(T,2);
FDATABASE(id).(regionlow)=mean(LOW);
FDATABASE(id).(regionhigh)=mean(HIGH);
FDATABASE(id).(regionamp)=mean(AMP);
FDATABASE(id).(regionframeij)=FRAMEij;
FDATABASE(id).(regiondur)=mean(DURij(:,1));
FDATABASE(id).(regionlowmatrix)=RLOWROI;
FDATABASE(id).(regionhighmatrix)=RHIGHROI;
FDATABASE(id).(regiondurmatrix)=DUROI;
k=round(mean(KS));
%FDATABASE(id).STIMTRACE=ST(k-round(totalscans/2)+1:k+round(totalscans/2));


% save database
if size(FDATABASE,2)>1
    PREVLOW=[];PREVHIGH=[];TR=[];
    for i=1:size(FDATABASE,2)-1
        PREVLOW=[PREVLOW;FDATABASE(i).(regionlow)];
        PREVHIGH=[PREVHIGH;FDATABASE(i).(regionhigh)];
        TR=[TR,FDATABASE(i).(regionmeantrace)];
    end 
    figure;plot(TR,'k');hold on;
    plot(FDATABASE(id).(regionmeantrace),'b')
    fprintf(['mean low: ',num2str(mean(PREVLOW)),' mean high:',num2str(mean(PREVHIGH)),'\n']);
    fprintf(['new low: ',num2str(FDATABASE(id).(regionlow)),' new high:',num2str(FDATABASE(id).(regionhigh)),'\n']);
    plot(1:length(TR),mean(PREVLOW),'r');plot(1:length(TR),mean(PREVHIGH),'r');
    plot(round(totalscans/2),mean(PREVLOW),'ro');plot(round(totalscans/2),mean(PREVHIGH),'ro');
    hold off;
end

% save to harddrive
FDATABASE(:).name
asksave=input('Save data <RETURN>');

if isempty(asksave)==1
    if exist('lastname')==0
        lastname=('DATABASE');
    end
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save DATABASE',lastname);
end
save([PathName,FileName],'FDATABASE');
FDATABASE(id)
