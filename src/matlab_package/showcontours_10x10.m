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
%range=70;%256
%dt=1.3;%256


%range=275;
%range=300;%range of times
dt=5;%isochronal separation (10 frames = 5 ms)
range=400;%New (ms; can change 300-400)

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

setframeumwidth=10;%frame unit length (40)[um]. If zero, use entire frame
setframeumheight=10;

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
%% select zoom area and plot new contours
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
            ax=FRAMECENTER(2);ay=FRAMECENTER(1);
        else
            fprintf('select center of frame \n');
            [ax,ay]=ginput(1);
            FRAMECENTER(1)=floor(ay);FRAMECENTER(2)=floor(ax);%pixel
            %define frame width and height   
            framewidth=round(setframeumwidth/pixelcalfactor_x);
            frameheight=round(setframeumheight/pixelcalfactor_y);
            FrameXLim=[round(ax-framewidth/2),round(ax-framewidth/2)+framewidth];
            FrameYLim=[round(ay-frameheight/2),round(ay-frameheight/2)+frameheight];
            
        end
    end
 
 
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
%% define zoom contour plot
caxis(COLORRANGE);
[CONTZOOM,hzoom]=contour(isoczoom,INTERPX,INTERPY,CONTM,[lowlimit:dt:highlimit],'LineColor','k');
set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
%set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
set(hzoom,'Fill',fillcontour,'LineWidth',lineweightzoom);
%% apply colorrange
axes(isoc)
colormap(jet(numcolors));
caxis(COLORRANGE);

axes(isoczoom)
colormap(jet(numcolors));
caxis(COLORRANGE);
%% clean contour plots
contourlines=get(h,'Children');
contourlineszoom=get(hzoom,'Children');
defaultminlength=100;
askminlength=defaultminlength;

while isempty(askminlength)==0
    minlength=askminlength;
    offcontour=[];oncontour=[];
    for i=1:length(contourlines)
        cx=get(contourlines(i),'XData');
        cy=get(contourlines(i),'YData');
        if size(cx,1)<=minlength
            offcontour=[offcontour;i];
        else
            oncontour=[oncontour;i];
        end
        
    end
    % toggle visibility of selected contour lines
    set(contourlines(offcontour),'visible','off');
    set(contourlineszoom(offcontour),'visible','off');
    set(contourlines(oncontour),'visible','on');
    set(contourlineszoom(oncontour),'visible','on');
    %ask to change minlength
    %askminlength=input(['minimum contour length [',num2str(minlength),']>']);
    askminlength=[];
end
%% FIT PIXELS: INTERPOLATE SELECTED AREA
if newfit==1 && vvectors==1
%% Estimate time interval for fit    
    if exist('D')==0 || (newestvel==1 && isempty(D)==0) || isempty(D)==1
        fprintf('click isochrones to estimate time window. ESC to exit.\n')
        axes(isoczoom)
        D=[];
        button=0;
        while button~=27
            if button~=27
                [AX,AY,p]=ginput(2);
                button=p(1);
                    if button ~=27
                        D=([D,[AX(2)-AX(1);AY(2)-AY(1)]])
                    end
            end
        end
        %calculate propagation velocities from difference vectors
        DISTANCE=[];
        for k=1:size(D,2)
            DISTANCE=[DISTANCE;norm(D(:,k))];
        end
        %estimate conduction velocity for fot
        estvel=mean(DISTANCE)/dt;%[pixels/frame]
        v=estvel*pixelcalfactor_x*1.0e-3*scanrate;%[mm/s]        
    end
    fprintf(['mean velocity [mm/s]:',num2str(v),'\n']);
%% restrict ROI within FRAME
if cv==1 && showframe==1
    askpolygon=input('Define closed polygon area in ROI? (1=YES)');
    if askpolygon==1
        axes(isoczoom);
        polygonhandle=impoly(isoczoom);
        POLYGONPOSITIONS=wait(polygonhandle);
        %POLYGONPOSITIONS are given relative to entire camera frame
    end
end
%% fit function
%% open Matlab Pool
%parallel processing tools
% p=gcp
% if isempty(p)==1
%     p=gcp;
%     poolsize = p.NumWorkers;
%     parpool(poolsize);
% end
%% CV measurement
[VFRAMEi,VFRAMEj,RMSEFRAME,XFRAMEINTERP_0,YFRAMEINTERP_0,SIGNALFRAMEINTERP_0]=fitpixels3(C,FrameXLim,FrameYLim,estvel,SIGNAL);
%VFRAMEi: conduction velocity in i-direction (line) [pixel/frame]
%VFRAMEj: conduction velocity in j-direction (column) [pixel/frame]
%XFRAMEINTERP_0: X coordinates of frame used for CV calculation
%YFRAMEINTERP_0: Y coordinates of frame
%% Litte fit statistics
%framepixels: total number of active pixels in frame
framepixels=0;
rmsepixels=0;
RMSELIST=[];%to determine smallest and largest RMSE
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            framepixels=framepixels+1;
            if 0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse
                rmsepixels=rmsepixels+1;
                RMSELIST=[RMSELIST;RMSEFRAME(i,j)];
            end
        end
    end
end
SIGNALFRAME=SIGNAL(FrameYLim(1):FrameYLim(end),FrameXLim(1):FrameXLim(end));%original SIGNALS
framearea=floor(bwarea(SIGNALFRAME))*pixelcalfactor_x*pixelcalfactor_y;%estimate of tissue area in frame
area=floor(bwarea(SIGNAL))*pixelcalfactor_x*pixelcalfactor_y;%estimate total excitable tissue area
fprintf(['number of active pixels available for fit: ',num2str(framepixels),'\n']);
fprintf(['active pixels with RMSE<',num2str(maxrmse),': ',num2str(rmsepixels),'\n']);
fprintf(['max RMSE: ',num2str(max(RMSELIST)),'\n']);
fprintf(['min RMSE: ',num2str(min(RMSELIST)),'\n']);
%% determine ROI polygon coordinates in interpolated frame space
%make sure that POLYROI from last run is wiped out
clear POLYROI
if askpolygon==1
    %1. determine interpolation factor used in fitfunction:
    %size of the interpolated frame/size of the original frame
    originalframex=FrameXLim(2)-FrameXLim(1)+1;
    interpolatedframex=size(SIGNALFRAMEINTERP_0,2);
    pixelinterval=1/(round(interpolatedframex/originalframex));
    %2. calculate binary image of the polygon in the interpolated
    %frame space
    POLYROI=zeros(size(SIGNALFRAMEINTERP_0));
    X=(POLYGONPOSITIONS(:,1)-FrameXLim(1)+1)/pixelinterval;
    Y=(POLYGONPOSITIONS(:,2)-FrameYLim(1)+1)/pixelinterval;
    M=size(POLYROI,1);N=size(POLYROI,2);
    POLYROI=poly2mask(X,Y,M,N);
end
%% calibrate velocities
%DEFINITIONS: VELVECTOR=[VFRAMEj(i,j);VFRAMEi(i,j)];
CALVFRAMEi=VFRAMEi*scanrate*pixelcalfactor_y*1.0e-3;%[mm/s]
CALVFRAMEj=VFRAMEj*scanrate*pixelcalfactor_x*1.0e-3;%[mm/s]
%% filter vectors
%determine mean and standard deviation of all velocities with RMSE<maxrmse
VLIST=[];
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1 && (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
            VLIST=[VLIST;norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2)];
        end
    end
end
meanvelall=mean(VLIST);
stdvelall=std(VLIST);
%filter vectors
VELINTERVAL=[meanvelall-5*stdvelall,meanvelall+5*stdvelall];
VPLOTBIN=zeros(size(SIGNALFRAMEINTERP_0));%binary matrix of filtered pixels
VPLOTLIST=[];%list of filtered velocities
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            if (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
                v=norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2);
                %if VELINTERVAL(1)<v && v<VELINTERVAL(2)
                if v>0
                    VPLOTBIN(i,j)=1;
                    VPLOTLIST=[VPLOTLIST;v];
                end
            end
        end
    end
end
                    

%show numerical results
meanvelocity=mean(VPLOTLIST);
stdvelocity=std(VPLOTLIST);

fprintf(['mean velocity of selected vectors: (',num2str(meanvelocity),' +- ',num2str(stdvelocity/sqrt(length(VPLOTLIST))),') mm/s \n']);
fprintf(['number of selected vectors: ',num2str(length(VPLOTLIST)),' out of total number: ',num2str(length(VLIST)),'\n']);
fprintf(['window size (w x h) :  ',num2str(frameumwidth),' x ',num2str(frameumheight),' um \n']);
%% plot velocity vectors
%plotvectors=0;
if plotvectors==1
    figure(heartfigurezoom)
    hold on
    for i=1:1:size(VPLOTBIN,1)
        for j=1:1:size(VPLOTBIN,2)
            if VPLOTBIN(i,j)==1
                x=XFRAMEINTERP_0(1,j);
                y=YFRAMEINTERP_0(i,1);
                vx=CALVFRAMEj(i,j);
                vy=CALVFRAMEi(i,j);
                vxn=vx/norm([vx,vy],2);
                vyn=vy/norm([vx,vy],2);
                if normvectorscale==0
                    %plot vector [vx,vy] at position [x,y] NOT NORMALIZED
                    quiver(x,y,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
                else
                    %plot vector [vx,vy] at position [x,y] NORMALIZED
                    quiver(x,y,vxn,vyn,normvectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
                end
            end
        end
    end
    hold off
end
%% angle unit vector   
DEFAULTPOINTS=[ACTIVATIONLIST(first,2),ACTIVATIONLIST(first,1);ACTIVATIONLIST(last,2),ACTIVATIONLIST(last,1)];
%DEFAULTPOINTS(2,:)=[46,24];
if exist('DV','var')==0
    %Default vector
    DV=(DEFAULTPOINTS(2,:)-DEFAULTPOINTS(1,:))';
end

if angleunitvectorinput==1
    %display default points
    axes(isoc)
    if exist('INPUTPOINTS','var')==1
        rectangle('Position',[INPUTPOINTS(1,1),INPUTPOINTS(1,2),1,1],'FaceColor','m');
        rectangle('Position',[INPUTPOINTS(2,1),INPUTPOINTS(2,2),1,1],'FaceColor','y'); 
    else
        rectangle('Position',[DEFAULTPOINTS(1,1),DEFAULTPOINTS(1,2),1,1],'FaceColor','b','EdgeColor','w');  
        rectangle('Position',[DEFAULTPOINTS(2,1),DEFAULTPOINTS(2,2),1,1],'FaceColor','r','EdgeColor','w');  
    end
    keepvector=input('Keep propagation direction [RETURN]>');
    if isempty(keepvector)==0
        %define new direction of propagation
        fprintf('select propagation direction (2 clicks)');
        axes(isoc);
        INPUTPOINTS=floor(ginput(2));
        %plot new points
        rectangle('Position',[INPUTPOINTS(1,1),INPUTPOINTS(1,2),1,1],'FaceColor','m');
        rectangle('Position',[INPUTPOINTS(2,1),INPUTPOINTS(2,2),1,1],'FaceColor','y'); 
        %calculate normalized difference vector
        DV=(INPUTPOINTS(2,:)-INPUTPOINTS(1,:))';
    end
end 
%angle unit vector
angleunitvector=DV/norm(DV,2);
%% plot mean vector
if plotsumvector==1;%plot sum velocity vector
    %calculate mean velocity vector
    VECTORLISTxy=[];%to calculate mean vector
    for i=1:size(VPLOTBIN,1)
        for j=1:size(VPLOTBIN,2)
            if VPLOTBIN(i,j)==1
                VECTORLISTxy=[VECTORLISTxy;[CALVFRAMEj(i,j),CALVFRAMEi(i,j)]];
            end
        end
    end
    MEANVECTOR=[sum(VECTORLISTxy(:,1)),sum(VECTORLISTxy(:,2))]./length(VECTORLISTxy);
              
    %plot sum vector into main isochronal map
    figure(heartfigure);hold on
    R=[FRAMECENTER(2)-floor(MEANVECTOR(1)/2),FRAMECENTER(1)-floor(MEANVECTOR(2)/2)];
    quiver(R(1),R(2),MEANVECTOR(1),MEANVECTOR(2),sumvectorscale,'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
    
    if angleunitvectorinput>0
        %plot angle unit vector into main isochronal map
        quiver(R(1),R(2),angleunitvector(1)*norm(MEANVECTOR,2),angleunitvector(2)*norm(MEANVECTOR,2),sumvectorscale,'LineWidth',sumvectorwidth,'Color',[0.5,0.5,0.5],'MaxHeadSize',1);
    end
        
    %plot sum vector into frame
    %figure(heartfigurezoom);hold on
    %quiver(FRAMECENTER(2)-floor(SUMVECTORZOOM(2)/2),FRAMECENTER(1)-floor(SUMVECTORZOOM(1)/2),SUMVECTORZOOM(2),SUMVECTORZOOM(1),0,...
    %    'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
    %hold off
    %plot example vector to confirm correct angle calculation
    
    if exist('PLOTVECTOR','var')==1
        %search for index in HEART SPACE
        figure(heartfigurezoom)
        hold on
        [ax,ay]=ginput(1);PLOTVECTOR=[ay,ax];
        i=1;while YFRAMEINTERP_0(i,1)<PLOTVECTOR(1);i=i+1;end
        j=1;while XFRAMEINTERP_0(1,j)<PLOTVECTOR(2);j=j+1;end
        v0=[VFRAMEj(i,j)*pixelcalfactor_j;VFRAMEi(i,j)*pixelcalfactor_i]*scanrate*1.0e-3*vectorscale;
        %plot vector into zoom map
        quiver(XFRAMEINTERP_0(1,j),YFRAMEINTERP_0(i,1),v0(1),v0(2),'LineWidth',sumvectorwidth,'Color',[1 0 0],'MaxHeadSize',1);
        %plot angle unit vector
        quiver(XFRAMEINTERP_0(1,j),YFRAMEINTERP_0(i,1),angleunitvector(1)*norm(v0,2),angleunitvector(2)*norm(v0,2),0,'LineWidth',sumvectorwidth,'Color',[0.5,0.5,0.5],'MaxHeadSize',1);
        fprintf(['angle:',num2str(ANGLES(i,j)),' degrees.\n']);
    end
     
end

%save sum vector
if exist('NSUMVECTOR')==1 && sframe>0
    asksavevel=input('save sum vector(1=yes)?');
    if isempty(asksavevel)==1 || asksavevel==1
        NSUMVECTOR=[NSUMVECTOR;SUMVECTOR];  
    else
        if asksaveframe==1 && size(NFRAME,1)>size(NSUMVECTOR,1)
            NFRAME(end,:)=[];%delete last frame
        end
    end
end
%% end of velocity calculation
end
%% add colorbars
%COLORBAR

caxis(COLORRANGE);
CLim=(COLORRANGE);

if showcolorbar==1
    axes(isoc);
    
    YTick=[COLORRANGE(1):division:COLORRANGE(2)];
    YTickLabel=[];
    for s=1:length(YTick)
        YTickLabel=[YTickLabel,{[num2str(round((YTick(s)-lowlimit)*1.0e3/scanrate)),' ms']}];
    end
    cbar=colorbar('position',colorbarposition,'YTick',YTick,'YTickLabel',YTickLabel,'Box','on','YColor','k','XColor','k','FontWeight','bold','YAxisLocation','right','Layer','top','CLim',CLim);
    set(cbar,'YTickMode','manual');    
    set(cbar,'YTick',YTick,'YTickLabel',YTickLabel)
end

axes(isoczoom)
colormap(jet(numcolors));
caxis(COLORRANGE);
CLim=(COLORRANGE);

% if showcolorbar==1
%     cbarzoom=colorbar('peer',isoczoom,'position',[0.85 0.02 0.05 0.3],'YTick',YTick,'YTickLabel',YTickLabel,'Box','on','YColor','w','XColor','w','FontWeight','bold','YAxisLocation','right','Layer','top','CLim',CLim);
%     set(cbarzoom,'YTickMode','manual');
%     %Textbox
%     if map>1
%         text(FrameXLim(1)+round(framewidth)-5,FrameYLim(1)+round(frameheight)-1,['\DeltaR=',num2str(dt),''],'color','w','FontWeight','bold');
%     else
%         text(FrameXLim(1)+round(framewidth)-5,FrameYLim(1)+round(frameheight)-1,['\Deltat=',num2str(dt/scanrate*1.0e3),' ms'],'color','w','FontWeight','bold');
%     end
% end
%% export images
askexport=input('export figures and velocities [RETURN]');

if exist('imagepath')==0 
    imagepath=datapath;
else
    if length(imagepath)<2
        imagepath=datapath;
    end
end
imagepath=uigetdir(imagepath,'SELECT FOLDER FOR IMAGES'); 
imagepath=[imagepath,'\'];
description=input('DESCRIPTION>','s');
basefile=[stackfile(1:5),'-',description,'-delay',num2str(lowdelay),'-range',num2str(range),'-WIEN',num2str(SPATWIEN(1)),'-MED',num2str(SPATMED(1))];

if isempty(askexport)==1
    
    %WRITE CONTOUR PLOT
    figure(heartfigure)
    %get current image and save as movie frame
    F=getframe(gcf);
    %convert movieframe back into image
    [IMAGE,imagemap]=frame2im(F);
    %crop image to remove black border
    IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
    %IMAGEC=IMAGE;
    %outputresolution [width,height]
    EXPORTIMAGE=imresize(IMAGEC,[outputresolution(2),outputresolution(1)],'bicubic');
    %write image to file
    p=[imagepath,basefile,'-CONTOUR.tif'];
    imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
       
    %WRITE ZOOM CONTOUR PLOT
    figure(heartfigurezoom)
    %get current image and save as movie frame
    F=getframe(gcf);
    %convert movieframe back into image
    [IMAGE,imagemap]=frame2im(F);
    %crop image to remove black border
    IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
    %IMAGEC=IMAGE;
    %output aspectratio of zoom image depends on frame size
    %outputresolutionzoom [width,height]
    EXPORTIMAGE=imresize(IMAGEC,[outputresolutionzoom(2),outputresolutionzoom(1)],'bicubic');
    %write image to file
    p=[imagepath,basefile,'-ZOOMumxy',num2str(round(frameumwidth)),'x',num2str(round(frameumheight)),'.tif'];
    imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);

    %WRITE HEART IMAGE
    if exist('RGBSHIFTB')==1
        axes(heartimageaxes)
        hold on
        if showimageframe>0
            rectangle('Position',fresizefactor*[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--');
        end
        set(heartimageaxes,'TickLength',[0 0]);
        %get current image and save as movie frame
        F=getframe(gcf);
        %convert movieframe back into image
        [IMAGE,imagemap]=frame2im(F);
        %crop image to remove black border
        IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
        %IMAGEC=IMAGE;
        %EXPORTIMAGE=imresize(IMAGEC,2,'bicubic');
        EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
        %write image to file
        p=[imagepath,basefile,'-HEART.tif'];
        imwrite(EXPORTIMAGE,p,'tif','Compression','none');
    end
end