%% user input
%colormap of divergence field  
%DATA MATRIX
C=cvdivergence(VFRAMEj,VFRAMEi);
SIGNAL = SIGNALFRAMEINTERP_0;

showframe=1;%1=show frame in isochronal map

%clear FRAMECENTER
framewidth=300;
frameheight=300;

opengl software;
SPATWIEN=[2 2];
SPATMED=[2 2];
delborder=3;%delete pixels at border to eliminate filtering artifacts
outputresolution=[1024,1024];%output image resolution [height length]

%color limits
highlimit=0.8;
lowlimit=-highlimit;

filename = [stackfile(1:end-4),'-DIV_limits',num2str(highlimit)];

resizefactor=10;

%frame
%NFRAME=[];%comment to save frame coordinates
sframe=0;
framecolor=[1 0 0];%color of frame
frameweight=2;
showcenter=0;%1=show only frame center
%% filter data matrix
if sum(SPATWIEN)==0 || isempty(SPATWIEN)==1
    WIENC=C;
else
    WIENC=wiener2(C,SPATWIEN);
end

if sum(SPATMED)==0 || isempty(SPATMED)==1
    MEDC=WIENC;
else
    MEDC=medfilt2(WIENC,SPATMED);
end
%% remove border pixels

BORDER=zeros(size(SIGNAL));
if delborder>0
    for j=1:delborder
        B=bwboundaries(SIGNAL,'noholes');
        BCOORDINATES=B{1};%get boundary coordinates
        %assemble boundary image
        BOUNDARY=zeros(size(SIGNAL));
        for i=1:length(BCOORDINATES)
            BOUNDARY(BCOORDINATES(i,1),BCOORDINATES(i,2))=1;
        end
        %add boundary to BORDER
        BORDER(BOUNDARY==1)=1;
        %shrink SIGNAL by new BOUNDARY
        SIGNAL(BOUNDARY==1)=0;
    end
end

%replace matrix by filtered matrix
D=MEDC;
%% Data range
%determine data range
DLIST=[];
for i=1:size(D,1)
    for j=1:size(D,2)
        if SIGNAL(i,j)==1 && D(i,j)>0
            DLIST=[DLIST;[i,j,D(i,j)]];
        end
    end
end
[B,IX]=sort(DLIST(:,3));
SORTDLIST=zeros(size(DLIST));
for i=1:length(IX)
    SORTDLIST(i,:)=DLIST(IX(i),:);
end
%% map data in [0 1] interval for color coding
MAPD=zeros(size(D));
%Normalize DATA matrix
for i=1:size(D,1)
    for j=1:size(D,2)
        if SIGNAL(i,j)>0
            if D(i,j)>highlimit
                MAPD(i,j)=1;
            elseif lowlimit<=D(i,j) && D(i,j)<=highlimit
                MAPD(i,j)=(D(i,j)-lowlimit)/(highlimit-lowlimit);
            else
                MAPD(i,j)=0;
            end
        end
    end
end
%% Resize data matrix
%resize SIGNAL matrix
RSIGNAL=imresize(SIGNAL,resizefactor,'bicubic');
RSIGNAL(RSIGNAL<0.5)=0;RSIGNAL(RSIGNAL>=0.5)=1;
%resize DATA matrix
RMAPD=imresize(MAPD,resizefactor,'bicubic');
RMAPD(RSIGNAL<1)=0;%resizing can lead to values outside the original 
RMAPD(RMAPD>1)=1;    
%convert to RGB image
[IMAPD,mapMAPD]=gray2ind(RMAPD,256);
RGBMAPD=ind2rgb(IMAPD,jet(256));

%set color for all pixels outside SIGNAL area white
for i=1:size(RGBMAPD,1)
    for j=1:size(RGBMAPD,2)
        if RSIGNAL(i,j)<1
            RGBMAPD(i,j,:)=[1,1,1];
        end
    end
end

%% create figure and axis

mapfig=figure('Name',filename,'MenuBar','none','Units','normalized',...
'Position',[0.25 0.25 0.5 0.5],'Color','w','Visible','off');
mapaxes=axes('Units','normalized','Position',[0 0 1 1],'Visible','on');

%draw images
axes(mapaxes)
image(RGBMAPD)
set(mapaxes,'TickLength',[0 0])
set(mapfig,'visible','on');
set(mapaxes,'YDIR','reverse');
set(mapaxes,'DataAspectRatio',[1 1 1]);
set(mapaxes,'Visible','off');

%% select zoom area and plot new colormap
%define width of frame
axes(mapaxes)

%select center of frame
if exist('FRAMECENTER')==1 && isempty(FRAMECENTER)==0
    ax=FRAMECENTER(2);ay=FRAMECENTER(1);
else
    fprintf('select center of frame \n');
    [ax,ay]=ginput(1);
    FRAMECENTER(1)=round(ay);FRAMECENTER(2)=round(ax);
end

FrameXLim=round([ax,ax+framewidth]-framewidth/2);
FrameYLim=round([ay,ay+frameheight]-frameheight/2);
        
%show frame
if showframe==1
    rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
end

if showcenter==1
    rectangle('Position',[FRAMECENTER(2),FRAMECENTER(1),1,1],'LineWidth',2,'EdgeColor',framecolor)
end
%% create zoomn figure and axis

mapfigzoom=figure('Name',filename,'MenuBar','none','Units','normalized',...
'Position',[0.25 0.25 0.5 0.5],'Color','w','Visible','off');
mapaxeszoom=axes('Units','normalized','Position',[0 0 1 1],'Visible','on');

%draw images
axes(mapaxeszoom)
image(RGBMAPD)

%set new limits
set(mapaxeszoom,'XLim',FrameXLim,'YLim',FrameYLim);


set(mapaxeszoom,'TickLength',[0 0])
set(mapfigzoom,'visible','on');
set(mapaxeszoom,'YDIR','reverse');
set(mapaxeszoom,'DataAspectRatio',[1 1 1]);
set(mapaxeszoom,'Visible','off');

%% export images
%construct file names
if exist('imagepath')==0 || isempty(imagepath)==1 || length(imagepath)<3
    imagepath=datapath;
end
imagepath=uigetdir(imagepath,'pick path to save images');imagepath=[imagepath,'\'];
description=input('DESCRIPTION>','s');

filemap=[filename,'-',description,'.tif'];
filemapzoom=[filename,'-',description,'_Z.tif'];

%get current image and save as movie frame
%MAP
imagefile=[imagepath,filemap];
imagefilezoom=[imagepath,filemapzoom];

axes(mapaxes)
set(mapaxes,'DataAspectratioMode','auto')
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%crop image to remove black border
IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
%IMAGEC=IMAGE;
%original aspect ratio
aspectr=size(RGBMAPD,1)/size(RGBMAPD,2);
EXPORTIMAGE=imresize(IMAGEC,[outputresolution(1)*aspectr,outputresolution(2)],'bicubic');
%write image to file
imwrite(EXPORTIMAGE,imagefile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
set(mapaxes,'DataASpectratioMode','manual');

axes(mapaxeszoom)
set(mapaxeszoom,'DataAspectratioMode','auto');
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%crop image to remove black border
IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
%IMAGEC=IMAGE;
aspectr=frameheight/framewidth;
EXPORTIMAGE=imresize(IMAGEC,[outputresolution(1)*aspectr,outputresolution(2)],'bicubic');
%write image to file
imwrite(EXPORTIMAGE,imagefilezoom,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
set(mapaxes,'DataASpectratioMode','manual');