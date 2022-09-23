% overlay data matrix with background image
%% user input
DATA=MED1APDMATRIX;%DATA MATRIX

SIGNALS=ASIGNALPIXELS;SIGNALS(SCATTER==0)=0;
%outline of SIGNALS image
BOUNDARYPIXELS=bwboundaries(SIGNALS);
B=BOUNDARYPIXELS{1};
BOUNDARY=zeros(size(SIGNALS));
for i=1:length(B)
    BOUNDARY(B(i,1),B(i,2))=1;
end
%new SIGNALS matrix with boundary removed
SIGNALSB=SIGNALS;SIGNALSB(BOUNDARY>0)=0;

%determine data range
DATALIST=[];
for i=1:size(DATA,1)
    for j=1:size(DATA,2)
        if SIGNALS(i,j)==1 && DATA(i,j)>0
            DATALIST=[DATALIST;[i,j,DATA(i,j)]];
        end
    end
end
[B,IX]=sort(DATALIST(:,3));
SORTDATALIST=zeros(size(DATALIST));
for i=1:length(IX)
    SORTDATALIST(i,:)=DATALIST(IX(i),:);
end


%set range
RANGE=[50 650];%map RANGE INTERVAL into [0 1];
mindisplay=100;
displayscalebars=0;
transparency=0.7;
B=SHIFTB;%background image
resizefactor=fresizefactor;
%% prepare data frame
MAPDATA=zeros(size(DATA));

%Normalize DATA matrix

for i=1:size(DATA,1)
    for j=1:size(DATA,2)
        if SIGNALS(i,j)>0
            if DATA(i,j)<=min(RANGE)
                MAPDATA(i,j)=0;
            end
            if DATA(i,j)>=max(RANGE)
                MAPDATA(i,j)=1;
            end
            if min(RANGE)<DATA(i,j) && DATA(i,j)<max(RANGE)
                MAPDATA(i,j)=(DATA(i,j)-min(RANGE))/(max(RANGE)-min(RANGE));
            end
        end
    end
end
%% prepare RGB images
%background image
BADJUST=imadjust(B);
[IB,mapB]=gray2ind(SHIFTB,256);%convert background to indexed image
RGBIMAGE=ind2rgb(IB,gray(256));
%signals matrix
RESSIGNALSB=imresize(SIGNALSB,resizefactor,'bicubic');
RESSIGNALSB(RESSIGNALSB<0.9)=0;RESSIGNALSB(RESSIGNALSB>=0.9)=1;
%data matrix
RESMAPDATA=imresize(MAPDATA,resizefactor,'bicubic');
%filter data matrix
RESMAPDATA(RESSIGNALSB<1)=0;%apply RESSIGNALSB matrix
RESMAPDATA(RESMAPDATA>1)=1;
[IMAPDATA,mapMAPDATA]=gray2ind(RESMAPDATA,2^16);
RGBMAPDATA=ind2rgb(IMAPDATA,jet(2^16));
%alphamatrix
minlevel=(mindisplay-min(RANGE))/(max(RANGE)-min(RANGE));%minimum colorlevel to show in image
ALPHA=double(ones(size(RESSIGNALSB))*transparency);
ALPHA(RESMAPDATA<=minlevel)=0;
%% create figure and axis
overlayfig=figure('Name','OVERLAY IMAGE','MenuBar','none','Units','normalized',...
    'Position',[0.25 0.25 0.5 0.55],'Color','w','Visible','off');
overlayaxes=axes('Units','normalized','Position',[0 0 1 1],'Visible','on','Drawmode','fast');

%draw images
image(RGBIMAGE);hold on
image(RGBMAPDATA,'alphadata',ALPHA);
set(overlayaxes,'TickLength',[0 0])
set(overlayfig,'visible','on');

%draw scalebars at the bottom right corner of image
scalefactor=160/72;%[um/pixels, calibration from 7/26/2007]
unitlength=20;%length of scalebar in um
scalebarlength=unitlength/scalefactor;%length of scalebar in pixels (80x80 image)
fractionscalebar=0.10;
unitstring=[num2str(unitlength),' \mum'];
scalebarxdata=resizefactor*(80-scalebarlength+1):resizefactor*fractionscalebar:resizefactor*80;
scalebarydata=resizefactor*80*ones(1,length(scalebarxdata));
xoffset=-10;
yoffset=-20;
if displayscalebars==1
    line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','w','LineWidth',2);
    line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','w','LineWidth',2);
end
    
%add unit length to horizontal scalebar
%text('Position',[resizefactor*70+xoffset+20,resizefactor*78+yoffset],'String',unitstring,'color','w','FontWeight','bold','FontSize',10)
%% additions to figure
%mark pixel positions previously defined in pixelmonitor
axes(overlayaxes)
if exist('N')==1
    TRACE=zeros(size(NORMDATA,3),size(N,1));
    for i=1:size(N,1)
        rectangle('Position',[(N(i,2)-0.5)*resizefactor,(N(i,1)-0.5)*resizefactor,resizefactor,resizefactor],'LineWidth',1.5,'EdgeColor','w');
    %save traces
    TRACE(:,i)=squeeze(NORMDATA(N(i,1),N(i,2),:));
    end
end
%% export image
%get current image and save as movie frame
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%crop image to remove black border
IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
%IMAGEC=IMAGE;
EXPORTIMAGE=imresize(IMAGEC,[768 1024],'bicubic');
%write image to file
[imagefile,imagepath]=uiputfile('*.tif','pick path');
description=input('DESCRIPTION>','s');
imagefile=[normfile(1:end-6),'-',description,'.tif'];
imagepathfile=[imagepath,imagefile];
imwrite(EXPORTIMAGE,imagepathfile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);