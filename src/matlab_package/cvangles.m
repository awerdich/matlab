function [ANGLES]=cvangles(angleunitvector,VFRAMEi,VFRAMEj,SIGNALFRAMEINTERP_0,imagepath,description)
%calculates divergence of conduction velocity field
%X=XFRAMEINTERP_0*pixelcalfactor*1.0e-3;% [mm]
%Y=YFRAMEINTERP_0*pixelcalfactor*1.0e-3;% [mm]
%VX=VFRAMEj*scanrate*pixelcalfactor*1.0e-3; %[pixel/frame*frame/s*um/pixel*1.0e-3mm/um]=[mm/s]
%VY=VFRAMEi*scanrate*pixelcalfactor*1.0e-3; %[mm/s]
%plot parameter
%[angleunitvector]=[x,y], image space!!!
range=180;
lowlimit=-range;
highlimit=range;
resizefactor=10;%resize divergence image for plotting
outputresolution=[1024,1024];%output image resolution [height length]
%% calculate angles relative to activation direction
%ALL CALCULATIONS IN IMAGE SPACE (x,y)!
%LIST=[LIST;[i,j,POL(i,j)]];
n=angleunitvector/norm(angleunitvector,2);
e=[1;0];%unit vector in x-direction
beta=acos(dot(e,n))*180/pi;
%convention: positive angles are measured counterclockwise
if n(2)>0 
    beta=-beta;
end
ANGLES=zeros(size(VFRAMEi));
ANGLESLIST=[];
for i=1:size(ANGLES,1)
    for j=1:size(ANGLES,2)
        if VFRAMEi(i,j)~=0 && VFRAMEj(i,j)~=0
            %normalized velocity vector in image space (x,y)
            v=[VFRAMEj(i,j),VFRAMEi(i,j)];
            v=v/norm(v,2);
            alpha=acos(dot(e,v))*180/pi;
            if v(2)>0
                alpha=-alpha;
            end
            %calculate the difference angle between alpha and n
            gamma=alpha-beta;
            %make sure that all angles are [-180,180]
            if gamma>180
                gamma=gamma-360;
            end
            if gamma<-180
                gamma=360+gamma;
            end
            ANGLES(i,j)=gamma;
            ANGLESLIST=[ANGLESLIST;gamma];
        else
            ANGLES(i,j)=-190;%not calculated angles
        end
    end
end
%% prepare to plot angles field if there is an angle field
if exist('SIGNALFRAMEINTERP_0','var')==1 && isempty(SIGNALFRAMEINTERP_0)==0
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
    %imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
end