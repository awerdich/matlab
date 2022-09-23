function [CVMAP]=cvcolormap(VFRAMEi,VFRAMEj,pixelcalfactor_i,pixelcalfactor_j,SIGNALFRAMEINTERP_0,imagepath,description);
%generate colormap of conduction velocities in ROI
%X=XFRAMEINTERP_0*pixelcalfactor*1.0e-3;% [mm]
%Y=YFRAMEINTERP_0*pixelcalfactor*1.0e-3;% [mm]
%VX=VFRAMEj*scanrate*pixelcalfactor*1.0e-3; %[pixel/frame*frame/s*um/pixel*1.0e-3mm/um]=[mm/s]
%VY=VFRAMEi*scanrate*pixelcalfactor*1.0e-3; %[mm/s]
%% plot parameter
lowlimit=0;
highlimit=10;
resizefactor=10;%resize divergence image for plotting
%outputresolution=resizefactor*[size(VFRAMEi,1),size(VFRAMEi,2)];%output image resolution [height length]
%% absolute conduction velocities in ROI
CVi=VFRAME


%% prepare to plot angles field
% Data range
ALIST=[];
for i=1:size(ANGLES,1)
    for j=1:size(ANGLES,2)
        if SIGNALFRAMEINTERP_0(i,j)==1 
            ALIST=[ALIST;[i,j,ANGLES(i,j)]];
        end
    end
end
SORTALIST=sortrows(ALIST,3);
%% map data in [0 1] interval for color coding
MAPD=zeros(size(ANGLES));
%Normalize DATA matrix
for i=1:size(ANGLES,1)
    for j=1:size(ANGLES,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            if ANGLES(i,j)>highlimit
                MAPD(i,j)=1;
            elseif lowlimit<=ANGLES(i,j) && ANGLES(i,j)<=highlimit
                MAPD(i,j)=(ANGLES(i,j)-lowlimit)/(highlimit-lowlimit);
            else
                MAPD(i,j)=0;
            end
        end
    end
end
%% Resize data matrix
%SIGNAL MATRIX
RSIGNAL=imresize(SIGNALFRAMEINTERP_0,resizefactor,'bicubic');
RSIGNAL(RSIGNAL<0.5)=0;RSIGNAL(RSIGNAL>=0.5)=1;
%data matrix
RMAPD=imresize(MAPD,resizefactor,'bicubic');
%filter data matrix
RMAPD(RSIGNAL==0)=0;%apply RESSIGNALSB matrix
RMAPD(RMAPD>1)=1;RMAPD(RMAPD<0)=0;%re-map values outside range    
[IMAPD,mapMAPD]=gray2ind(RMAPD,256);
RGBMAPD=ind2rgb(IMAPD,jet(256));
%color holes white
for i=1:size(RSIGNAL,1)
    for j=1:size(RSIGNAL,2)
        if RSIGNAL(i,j)==0
            RGBMAPD(i,j,:)=[1,1,1];
        end
    end
end
%% create figure and axis
P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
windowsize=0.75.*[maxwindowsize,maxwindowsize];%[width,height]*screensize
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/2);%distance of windows from bottom of screen
sep=round((screenwidth-2*windowsize(1)-2*offsetx)/2);%distance between windows


%create figure
mapfig = figure('Name','ANGLES','MenuBar','none','Units','pixels','Position',[offsetx,offsety,windowsize],'Color','w','Visible','on');

%set axes for mapfigure
figure(mapfig);
mapaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);


%draw images
axes(mapaxes)
image(RGBMAPD);

set(mapaxes,'TickLength',[0 0])
set(mapfig,'visible','on');
set(mapaxes,'YDIR','reverse');
set(mapaxes,'DataAspectRatio',[1 1 1]);
%% export image
basefile=['ANGLES-',description,'-RANGE ',num2str(range)];


    %WRITE CONTOUR PLOT
    figure(mapfig)
    %get current image and save as movie frame
    F=getframe(mapaxes);
    %convert movieframe back into image
    [IMAGE,imagemap]=frame2im(F);
    %crop image to remove black border
    IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
    %IMAGEC=IMAGE;
    EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
    %write image to file
    p=[imagepath,basefile,'.tif'];
    imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);