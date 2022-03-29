%% SETTINGS
%data input matrix
C=POL;
%save choice of input matrix
if min(all(C==POL))==1
    %C=POL
    pol1repol2=1;
else
    %C=REPOL
    pol1repol2=2;
end

clear FRAMECENTER; 
%close all
if exist('magnification','var')==0 || isempty(magnification)==1; mapconfig; end
%define window size based on actual screen size
P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
windowsize=0.75.*[maxwindowsize,maxwindowsize];%[width,height]*screensize
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/4);%distance of windows from bottom of screen
sep=round((screenwidth-2*windowsize(1)-2*offsetx)/2);%distance between windows

lineweight=2.5;%line width of contour plot
lineweightzoom=3.0;%line width of zoom contour plot
markframecenter=0;%mark pixel in frame;
showcolorbar=0;
showscalebar=0;%0 off 1 xonly 2 xy


lowdelay=-1;%time delay to plot after smallest value4
range=70;
dt=2;

%range=275;
%range=300;%range of times
%dt=10;%isochronal separation

matrixaspectratio=size(C,2)/size(C,1);%aspectratio of data input matrix [width/height]
outputresolution=[1024,round(1024/matrixaspectratio)];%output image resolution [width, height]
delborder=2;%number of border pixels to remove from map (FILTER artifacts)

%COLORBAR
lright=[0.75 0.02 0.05 0.3];%right side colorbar 'position'
lleft=[0.01 0.02 0.05 0.3];%left side colorbar 'position'
uleft=[0.01 0.68 0.05 0.3];
uright=[0.72 0.68 0.05 0.3];
colorbarposition=uleft;
numcolors=225;
colormap(jet(numcolors));

%FRAME
%NFRAME=[];%comment to save frame coordinates
sframe=1;%save frames and sum vectors

%setframeumwidth=25;%frame unit length (40)[um]. If zero, use entire frame
%setframeumheight=25;

%setframeumwidth=0;

setframeumwidth=30;%frame unit length (40)[um]. If zero, use entire frame
setframeumheight=30;

%setframeumwidth=1500;%frame unit length (40)[um]. If zero, use entire frame
%setframeumheight=1500;

%VECTORS
%NSUMVECTOR=[];%comment to save sum vectors
cv=1;%1:calculate CVs and show frame 0:show only isochronal map
%angle unit vector
if cv==1
    showframe=1;%1=show frame in isochronal map
    showimageframe=0;%1=show image frame
    vvectors=1;%1=show velocity vectors on contour plot
    plotvectors=1;%1=plot velocity vectors
    plotsumvector=0;%plot sum velocity vector
    angleunitvectorinput=1;%manually input angle unit vecotor 
    newfit=1;newestvel=1;
    fillcontour='on';
else
    showframe=0;%1=show frame in isochronal map
    showimageframe=0;%1=show image frame
    vvectors=0;%1=show velocity vectors on contour plot
    plotvectors=0;%1=plot velocity vectors
    plotsumvector=0;%plot sum velocity vector
    angleunitvectorinput=0;%manually input angle unit vecotor 
    newfit=0;newestvel=0;
    fillcontour='on';
end

% %perform new fit only if new area selected and if vectors are displayed
% if vvectors==1 && (exist('FRAMECENTER')==0 || isempty(FRAMECENTER)==1)
%     newfit=1;newestvel=1;
% end
%maxrmse=0.75;%maximum root-mean-square error allowed [frames]
maxrmse=4;
%EMBRYOS
vectorscale=0.2;%scaling factor of velocity vector field (usually 0.2)
sumvectorscale=2;
normvectorscale=0;%scaling of normalized vector 0:not normalized

%ADULT
%vectorscale=0.02;%scaling factor of velocity vector field (usually 0.2)
%sumvectorscale=0.25;

vectorwidth=1.5;%[pts] vector line width [2.5]
sumvectorwidth=4;
unitvectorwidth=4;
normsumvector=0;%normalizes sum vector so that all sumvectors have the same length

%PIXEL CALIBRATION FACTOR
unitlength=20;%length of scalebar in um

%COLORS
framecolor=[0 0 0];%color of frame
frameweight=3;
vectorcolor=[0 0 0];
sumvectorcolor=[1 0 0];
unitvectorcolor=[0.5 0.5 0.5];
%% process data matrix
SPATWIEN=[2 2];
SPATMED=[2 2];
%CONTOUR PLOT
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

%SIGNAL MATRIX
SIGNAL=ASIGNALPIXELS;
SIGNAL(SCATTER==0)=0;
SIGNAL(POL==0)=0;

%Calculate 1 perimeter
B1=bwperim(SIGNAL);
%remove one border pixel from SIGNAL
SIGNAL(B1==1)=0;

%replace C by filtered matrix
C2=C;%make a backup of original matrix
C=MEDC;C(SIGNAL==0)=0;

%replace second border by unfiltered data
B2=bwperim(SIGNAL);
for i=1:size(C,1)
    for j=1:size(C2,2)
        if B2(i,j)==1
            C(i,j)=C2(i,j);
        end
    end
end
%% delete signals outside the tissue area caused by light scattering
%number of ASIGNALPIXELS
numallsignal=0;
for i=1:size(SIGNAL,1)
    for j=1:size(SIGNAL,2)
        if SIGNAL(i,j)==1
            numallsignal=numallsignal+1;
        end
    end
end
%% find limits of times
TIMES=[];
for i=1:size(C,1)
    for j=1:size(C,2)
        if SIGNAL(i,j)>0
            TIMES=[TIMES;[C(i,j),i,j]];
        end
    end
end
%sort times
[S,SIDX]=sort(TIMES(:,1));
STIMES=TIMES(SIDX,:);
%% ISOCHRONE PARAMETERS
%ACTIVATION MAP
%alternatively, provide step size for colorbar
stepsize=10;
%TIME LIMITS
lowlimit=STIMES(1,1)+lowdelay;
highlimit=lowlimit+range;
COLORRANGE=[lowlimit,highlimit];
division=(highlimit-lowlimit)/2;
tickdivision=numcolors/(highlimit-lowlimit)*division;
%% create contourmatrix
%set range of all pixels between minframe and maxframe

for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j)<lowlimit
            C(i,j)=-1e8;
            %delete pixels with activation times below lowlimit from SIGNAL
            SIGNAL(i,j)=0;
        elseif C(i,j)>highlimit
            C(i,j)=1e8;
            SIGNAL(i,j)=0;
        end
    end
end
%% determine first and last pixel in ACTIVATIONLIST and define unitvector
k=1;while SIGNAL(ACTIVATIONLIST(k,1),ACTIVATIONLIST(k,2))==0;k=k+1;end;first=k+10;
k=size(ACTIVATIONLIST,1);while SIGNAL(ACTIVATIONLIST(k,1),ACTIVATIONLIST(k,2))==0;k=k-1;end;last=k-50;

%% prepare figures and axes
    
%define axes for picture display

heartfigure = figure('Name','CONTOUR PLOT','MenuBar','none','Units','pixels','Position',[offsetx,offsety,windowsize(1),windowsize(2)/matrixaspectratio],'Color','w','Visible','on');
heartfigurezoom = figure('Name','CONTOUR PLOT ZOOM','MenuBar','none','Units','pixels','Position',[offsetx+sep+windowsize(1),offsety,windowsize(1),windowsize(2)/matrixaspectratio],'Color','w','Visible','on');
heartimage = figure('Name','HEARTIMAGE','MenuBar','none','Units','pixels','Position',[offsetx+sep+windowsize(1),offsety,windowsize],'Color','w','Visible','on');

%set axes for figheart
figure(heartfigure);
isoc=axes('Position',[0 0 1 1],'Visible','on','YDir','reverse','TickLength',[0 0]);

figure(heartfigurezoom);
isoczoom=axes('Position',[0 0 1 1],'Visible','on','YDir','reverse','TickLength',[0 0]);

if exist('RGBSHIFTB')==1
    figure(heartimage)
    heartimageaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','TickLength',[0 0]);
    axes(heartimageaxes)
    image(RGBSHIFTB);
end

%% interpolate contour matrix
%prepare contour matrix
contresizefactor=10;
[X,Y]=meshgrid([1:size(C,2)],[1:size(C,1)]);%grid for 2D interpolation
[INTERPX,INTERPY]=meshgrid([1:1/contresizefactor:size(X,2)],[1:1/contresizefactor:size(Y,1)]);%interpolated grid
CONTM=interp2(X,Y,C,INTERPX,INTERPY,'linear');%interpolated contourmatrix

%% draw contour plot
axes(isoc);
caxis(COLORRANGE);
[CONT,h]=contour(isoc,INTERPX,INTERPY,CONTM,[lowlimit:dt:highlimit],'LineColor','k');
set(isoc,'YDir','reverse','Color','none','TickDir','out','XLim',[1,size(NORMDATA,2)],'YLim',[1,size(NORMDATA,1)]);
set(h,'Fill',fillcontour,'LineWidth',lineweight);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
%add scalebar
scalebarlength=unitlength/pixelcalfactor_x;%length of scalebar in pixels (80x80 image)
fractionscalebar=0.10;
unitstring=[num2str(unitlength),'\mum'];
scalebarxdata=(80.0-scalebarlength):fractionscalebar:80.0;
scalebarydata=80*ones(1,length(scalebarxdata));
%draw scalebars at the bottom right corner of image
xoffset=-1;
yoffset=-2;
if showscalebar>0
    line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','k','LineWidth',2);
    if showscalebar>1
        line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','k','LineWidth',2);
    end 
%add unit length to horizontal scalebar
text(70,78+yoffset,unitstring,'color','k','FontWeight','bold')
text(70,78+yoffset-3,['\Deltat=',num2str(dt/scanrate*1.0e3),'ms'],'color','k','FontWeight','bold');
end
%% Draw Contour Map
%select zoom area and plot new contours
%define width of frame
axes(isoc)
    if (setframeumwidth==0 || setframeumheight==0)
        if exist('FRAMECENTER','var')==0
            %use crop tool to define frame
            %convert plot into image
            F=getframe(isoc);
            [I,Map]=frame2im(F);
            displayfactor=2;%reduce to fit on screen
            RI=imresize(I,[size(Y,1),size(X,2)].*contresizefactor./displayfactor);
            %determine crop coordinates
            figure;
            [I2,RECT]=imcrop(RI);
            %convert rectangle back to pixel coordinates
            PRECT=round(RECT*displayfactor/contresizefactor);
            FRAMECENTER=[PRECT(2)+PRECT(4)/2,PRECT(1)+PRECT(3)/2];
        end
        framewidth=round(PRECT(3));
        frameheight=round(PRECT(4));
        setframeumwidth=framewidth*pixelcalfactor_x;
        setframeumheight=frameheight*pixelcalfactor_y;
        FrameXLim=[PRECT(1),PRECT(1)+framewidth];
        FrameYLim=[PRECT(2),PRECT(2)+frameheight];
    else
         %select center of frame
        if exist('FRAMECENTER')==1 && isempty(FRAMECENTER)==0
            ax(1)=FRAMECENTER(2);ay(1)=FRAMECENTER(1);
        else
            fprintf('select center of regions \n');
            [ax,ay]=ginput(8);
        end
            
    end
    [v,d] = size(ax);
    b = 0;
    while b < v
        b = b + 1;    
        FRAMECENTER(1)=floor(ay(b));FRAMECENTER(2)=floor(ax(b));%pixel
            %define frame width and height   
           
            framewidth=round(setframeumwidth/pixelcalfactor_x);
            frameheight=round(setframeumheight/pixelcalfactor_y);
            FrameXLim=[round(ax(b)-framewidth/2),round(ax(b)-framewidth/2)+framewidth];
            FrameYLim=[round(ay(b)-frameheight/2),round(ay(b)-frameheight/2)+frameheight];
            %end
 %check if Frame limits are consistent with array limit
 if FrameXLim(1)<1,framewidth=framewidth-(1-FrameXLim(1));FrameXLim(1)=1;end
 if FrameXLim(2)>size(C,2),framewidth=framewidth-(size(C,2)-FrameXLim(2));FrameXLim(2)=size(C,2);end
 if FrameYLim(1)<1,frameheight=frameheight-(1-FrameYLim(1));FrameYLim(1)=1;end
 if FrameYLim(2)>size(C,1),frameheight=frameheight-(size(C,1)-FrameYLim(2));FrameYLim(2)=size(C,1);end        
        
 %re-set final frame width and frame height using FrameXLim, FrameYLim
 framewidth=FrameXLim(2)-FrameXLim(1);
 frameheight=FrameYLim(2)-FrameYLim(1);
 frameumwidth=framewidth*pixelcalfactor_x; 
 frameumheight=frameheight*pixelcalfactor_y; 

 
 % add frame to mark selected area
axes(isoc);
if showframe==1
    set(isoc,'DrawMode','normal');
    rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
    %rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--');
end
    end
% save frame coordinates
if exist('NFRAME')>0 && sframe>0
    asksaveframe=input('save frame (1=yes)?');
    if asksaveframe==1
        NFRAME=[NFRAME;[FrameXLim,FrameYLim]];
   end
end 
% adjust aspectratio of zoom image
frameaspectratio=framewidth/frameheight;
WINDOWSIZEZOOM=zeros(size(windowsize));
if maxwindowsize<windowsize(2)/frameaspectratio
    WINDOWSIZEZOOM(2)=windowsize(2);
    WINDOWSIZEZOOM(1)=windowsize(2)*frameaspectratio;
else
    WINDOWSIZEZOOM(2)=windowsize(2)/frameaspectratio;
    WINDOWSIZEZOOM(1)=windowsize(2);
end
set(heartfigurezoom,'Position',[offsetx+sep+windowsize(1),offsety,WINDOWSIZEZOOM(1),WINDOWSIZEZOOM(2)]);

outputresolutionzoom=[outputresolution(1),round(outputresolution(2)/frameaspectratio)];%resolution of zoom image [width,height]
%%Select regions in a specific order with allowing to ommit some regions
%%Calculations for each region
