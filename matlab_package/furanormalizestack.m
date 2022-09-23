%normalize fura2 image stack for movie 
%% collect info
DATA=EPISODEFITR;
NORMDATA=zeros(size(DATA));
SIGNAL=ASIGNALPIXELS;
SIGNAL(SCATTER==0)=0;
fmax=0.4;%RED
fmin=0.0;%BLUE
%% Calculate Fura matrices
TMIN=zeros(size(DATA,1),size(DATA,2));
TMAX=zeros(size(TMIN));
AMPLIST=[];
for i=1:size(DATA,1)
    for j=1:size(DATA,2)
        if SIGNAL(i,j)==1
            PIXEL=squeeze(DATA(i,j,:));
            BASELINESCANS=squeeze(BSCANS(i,j,:));
            %BASELINESCANS=[1:20]';
            MAXIMUMSCANS=squeeze(MSCANS(i,j,:));
            TMIN(i,j)=mean(PIXEL(BASELINESCANS));
            TMAX(i,j)=mean(PIXEL(MAXIMUMSCANS));
            AMPLIST=[AMPLIST;TMAX(i,j)-TMIN(i,j)];
            %remove baseline offset
            OPIXEL=PIXEL-TMIN(i,j);
            NPIXEL=(OPIXEL-fmin)/(fmax-fmin);
            NORMDATA(i,j,:)=NPIXEL;
        end
    end
end