%% SETTINGS
close all
clear ACTIVATIONLIST
windowsize=[0.3,0.5];%[width,height]*screensize
%CONTOUR PLOT
C=POL;
first=50;%line number in ACTIVATIONLIST for the first active pixel coordinates
if exist('ACTIVATIONLIST')==1;last=size(ACTIVATIONLIST,1);end
lineweight=2;%line width of contour plot
lineweightzoom=2.4;%line width of zoom contour plot
markframecenter=0;%mark pixel in frame;
pixelinput=0;%input pixels and display in contour map
showcolorbar=0;
showscalebar=0;%0 off 1 xonly 2 xy
SPATWIEN=[3 3];
SPATMED=[1 1];
lowdelay=-9;%time delay to plot after smallest value
range=800;%range of times
dt=10;%isochronal separation
fillcontour='on';
outputresolution=[1024,1024];%output image resolution [height length]
delborder=2;%number of border pixels to remove from map (FILTER artifacts)

%COLORBAR
lright=[0.75 0.02 0.05 0.3];%right side colorbar 'position'
lleft=[0.01 0.02 0.05 0.3];%left side colorbar 'position'
uleft=[0.01 0.68 0.05 0.3];
uright=[0.72 0.68 0.05 0.3];
colorbarposition=uright;
numcolors=225;
colormap(jet(numcolors));

%FRAME
%NFRAME=[];%comment to save frame coordinates
sframe=1;%save frames and sum vectors

showframe=1;%1=show frame in isochronal map
showimageframe=0;%1=show image frame
framecenterline=first;
%FRAMECENTER=[ACTIVATIONLIST(framecenterline,1),ACTIVATIONLIST(framecenterline,2)];
%NFRAMECENTER=[];NSUMVECTOR=[];%save frame center and summary vectors for all subsequent measurements
setframeumwidth=35;%frame unit length (40)[um]. If zero, click corners of frame
setframeumheight=35;
% FRAMECENTER=DATABASE(1).a_windowcenter;
% setframeumwidth=round(DATABASE(2).av_windowx);
% setframeumheight=round(DATABASE(2).av_windowy);
% N=DATABASE(2).av_pixcoordinates;

%VECTORS
%NSUMVECTOR=[];%comment to save sum vectors
vvectors=1;%1=show velocity vectors on contour plot
plotvectors=1;%1=plot velocity vectors
plotsumvector=1;%plot sum velocity vector
newfit=0;newestvel=0;
%perform new fit only if new area selected and if vectors are displayed
if vvectors==1 && (exist('FRAMECENTER')==0 || isempty(FRAMECENTER)==1)
    newfit=1;newestvel=1;
end
maxrmse=0.7;%maximum root-mean-square error allowed [frames]
vectorscale=0.2;%scaling factor of velocity vector field v=0.075;av=1 and o=1
sumvectorscale=0.002;
vectorwidth=2;%[pts] vector line width
sumvectorwidth=3.5;
normsumvector=0;%normalizes sum vector so that all sumvectors have the same length

%PIXEL CALIBRATION FACTOR
pixelcalibration=16/7.201613;%REDSHIRT 80x80 [um/pixel] at 20x Objective lens, 0.5x sideport
magnification=20;%objective magnification
pixelcalfactor=pixelcalibration/magnification*20;
unitlength=20;%length of scalebar in um

%COLORS
framecolor=[1 0 0];%color of frame
frameweight=3.5;
vectorcolor=[0 0 0];
sumvectorcolor=[0 0 0];
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
%% replace boder of active pixels 
%combine ASIGNALPIXELS and SCATTER in SIGNAL
SIGNAL=ASIGNALPIXELS;
SIGNAL(SCATTER==0)=0;
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
%replace BORDER pixels with unfiltered data
for i=1:size(BORDER,1)
    for j=1:size(BORDER,2)
        if BORDER(i,j)==1
            MEDC(i,j)=C(i,j);
        end
    end
end
%replace contour matrix by filtered matrix
C=MEDC;
%restore SIGNAL matrix
% SIGNAL=ASIGNALPIXELS;
% SIGNAL(SCATTER==0)=0;
%% delete signals outside the tissue area caused by light scattering
%number of ASIGNALPIXELS
numallsignal=0;
for i=1:size(SIGNAL,1)
    for j=1:size(SIGNAL,2)
        if SIGNAL(i,j)==0
            C(i,j)=0;
        else
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
division=stepsize;
tickdivision=numcolors/(highlimit-lowlimit)*division;
%% create contourmatrix
%set range of all pixels between minframe and maxframe
M=C;%create copy of the original time matrix before making changes
C(SIGNAL==0)=0;
M(SIGNAL==0)=0;
for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j)<lowlimit
            C(i,j)=-1e8;
        elseif C(i,j)>highlimit
            C(i,j)=1e8;
        end
    end
end
%% prepare figures and axes
    
%define axes for picture display

heartfigure = figure('Name','CONTOUR PLOT','MenuBar','none','Units','normalized','Position',[0,0.4,windowsize],'Color','w','Visible','on');
heartfigurezoom = figure('Name','CONTOUR PLOT ZOOM','MenuBar','none','Units','normalized','Position',[0.5 0.4 windowsize],'Color','w','Visible','on');
heartimage = figure('Name','HEARTIMAGE','MenuBar','none','Units','normalized','Position',[0.5 0.4 windowsize],'Color','w','Visible','on');

%set axes for figheart
figure(heartfigure);
isoc=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);

figure(heartfigurezoom);
isoczoom=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);

if exist('RGBSHIFTB')==1
    figure(heartimage)
    heartimageaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','TickLength',[0 0]);
    axes(heartimageaxes)
    image(RGBSHIFTB);
end
%% interpolate contour matrix
%prepare contour matrix
gridsize=size(C,1);%pixelnumber, assuming a square array
[X,Y]=meshgrid([1:gridsize]);%grid for 2D interpolation
contresizefactor=10;
contfraction=(gridsize-1)/(contresizefactor*gridsize-1);%division of grid
[INTERPX,INTERPY]=meshgrid([1:contfraction:gridsize]);%interpolated grid
CONTM=interp2(X,Y,C,INTERPX,INTERPY,'linear');%interpolated contourmatrix
%% draw contour plot
axes(isoc);
caxis(COLORRANGE);
[CONT,h]=contour(isoc,INTERPX,INTERPY,CONTM,[lowlimit:dt:highlimit],'LineColor','k');
set(isoc,'YDir','reverse','Color','none','TickDir','out','XLim',[1,80],'YLim',[1,80]);
set(h,'Fill',fillcontour,'LineWidth',lineweight);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
%add scalebar
scalebarlength=unitlength/pixelcalfactor;%length of scalebar in pixels (80x80 image)
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
if map>1
    text(70,78+yoffset-3,['\DeltaR=',num2str(dt),''],'color','k','FontWeight','bold');
else
    text(70,78+yoffset-3,['\Deltat=',num2str(dt/scanrate*1.0e3),'ms'],'color','k','FontWeight','bold');
end
end
%%mark first and last activated pixels
if exist('ACTIVATIONLIST')==1
    rectangle('Position',[ACTIVATIONLIST(first,2)-0.5,ACTIVATIONLIST(first,1)-0.5,1,1],'FaceColor','g');  
    rectangle('Position',[ACTIVATIONLIST(last,2)-0.5,ACTIVATIONLIST(last,1)-0.5,1,1],'FaceColor','r');  
end
%% select zoom area and plot new contours
%define width of frame
axes(isoc)

    %select corners of frame
    if setframeumwidth==0 || setframeumheight==0
        fprintf('select two corners of frame \n');
        for k=1:2
            [ax(k),ay(k)]=ginput(1);
            %rectangle('Position',[ax(k)-0.5,ay(k)-0.5,1,1],'LineWidth',2,'EdgeColor','r')
        end    
        FrameXLim=round([min(ax),max(ax)]);FrameYLim=round([min(ay),max(ay)]);
        framewidth=FrameXLim(2)-FrameXLim(1);frameheight=FrameYLim(2)-FrameYLim(1);
        frameumwidth=round(framewidth*pixelcalfactor);frameumheight=round(frameheight*pixelcalfactor);
        FRAMECENTER=[mean(FrameYLim),mean(FrameXLim)];
    else
        %select center of frame
        if exist('FRAMECENTER')==1 && isempty(FRAMECENTER)==0
            ax=FRAMECENTER(2);ay=FRAMECENTER(1);
        else
            fprintf('select center of frame \n');
            [ax,ay]=ginput(1);
            FRAMECENTER(1)=ay;FRAMECENTER(2)=ax;
        end
        framewidth=setframeumwidth/pixelcalfactor;frameheight=setframeumheight/pixelcalfactor;    
        FrameXLim=round([ax,ax+framewidth]-framewidth/2);
        FrameYLim=round([ay,ay+frameheight]-frameheight/2);
        
        
        %check if Frame limits are consistent with array limit
        if FrameXLim(1)<1,framewidth=framewidth-(1-FrameXLim(1));FrameXLim(1)=1;end
        if FrameXLim(2)>size(C,2),framewidth=framewidth-(size(C,2)-FrameXLim(2));FrameXLim(2)=size(C,2);end
        if FrameYLim(1)<1,frameheight=frameheight-(1-FrameYLim(1));FrameYLim(1)=1;end
        if FrameYLim(2)>size(C,1),frameheight=frameheight-(size(C,1)-FrameYLim(2));FrameYLim(2)=size(C,1);end        
        
        frameumwidth=(FrameXLim(2)-FrameXLim(1))*pixelcalfactor; frameumheight=(FrameYLim(2)-FrameYLim(1))*pixelcalfactor;
    end
  
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
defaultminlength=200;
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
    % ask to change minlength
    askminlength=input(['minimum contour length [',num2str(minlength),']>']);
    %askminlength=[];
end
%% get pixels
if exist('pixelinput')==1 && pixelinput>0
    %pixel list in matrix coordinates (i,j)
    if exist('N')==0
        N=[];
        axes(isoczoom);
        fprintf(['mark points in map. <ENTER> when finished.']);
        [px,py]=ginput;
        for i=1:length(px)
            N(i,:)=round([py(i),px(i)]);
        end
    end
    %mark pixels in graph
    axes(isoczoom);
    for i=1:size(N,1)
        rectangle('Position',[N(i,2)-0.5,N(i,1)-0.5,1,1],'EdgeColor','r','LineWidth',1);
    end
    axes(isoc)
    for i=1:size(N,1)
        rectangle('Position',[N(i,2)-0.5,N(i,1)-0.5,1,1],'EdgeColor','r','LineWidth',1);
    end
end
%% FIT PIXELS: INTERPOLATE SELECTED AREA
if newfit==1
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
        v=estvel*pixelcalfactor*1.0e-3*scanrate;%[mm/s]        
    end
    fprintf(['mean velocity [mm/s]:',num2str(v),'\n']);

%% fit function
[VELOCITYLIST,framepixels]=fitpixels2(M,FrameXLim,FrameYLim,estvel,SIGNAL);
%VELVECTOR=pixelinterval*[Vi(i,j);Vj(i,j)];
VSTRING='[YFRAMEINTERP(i,1),XFRAMEINTERP(1,j),v,a,VELVECTOR(1),VELVECTOR(2),RMSE(i,j)] (matrix coordinates!!!)';        
% some fit statistics
%total number of pixels
fprintf(['number of active pixels used in fit: ',num2str(framepixels),'\n']);
fprintf(['successful fits:',num2str(size(VELOCITYLIST,1)),'\n']);
fprintf(['max RMSE: ',num2str(max(VELOCITYLIST(:,7))),'\n']);
fprintf(['min RMSE: ',num2str(min(VELOCITYLIST(:,7))),'\n']);
%% take only velocity vectors with RMSE<maxrmse
RVELOCITYLIST=[];%VELOCITIES with RMSE<maxrmse
for i=1:size(VELOCITYLIST,1)
    if VELOCITYLIST(i,7)<maxrmse
        RVELOCITYLIST=[RVELOCITYLIST;VELOCITYLIST(i,:)];
    end
end
fprintf(['successful fits with RMSE<',num2str(maxrmse),': ',num2str(size(RVELOCITYLIST,1)),'\n']);
fprintf(['fit success rate : ',num2str(size(RVELOCITYLIST,1)/framepixels),'\n']);
%% calibrate velocities
CALVELOCITYLIST=RVELOCITYLIST;
CALVELOCITYLIST(:,3)=RVELOCITYLIST(:,3)*scanrate*pixelcalfactor*1.0e-3;
CALVELOCITYLIST(:,5)=RVELOCITYLIST(:,5)*scanrate*pixelcalfactor*1.0e-3;
CALVELOCITYLIST(:,6)=RVELOCITYLIST(:,6)*scanrate*pixelcalfactor*1.0e-3;
end
if vvectors==1
%% select vectors
meanvelall=mean(CALVELOCITYLIST(:,3));
stdvelall=std(CALVELOCITYLIST(:,3));
meanangleall=mean(CALVELOCITYLIST(:,4));
stdangleall=std(CALVELOCITYLIST(:,4));
%filter vectors
VELINTERVAL=[0,meanvelall+2*stdvelall];
ANGINTERVAL=[0,360];
%ANGINTERVAL=[meanangleall-stdangleall,meanangleall+stdangleall];
%ANGINTERVAL=[meanangleall-2*stdangleall,meanangleall+2*stdangleall];
%ANGINTERVAL=[meanangleall-3*stdangleall,meanangleall+3*stdangleall];
VPLOT=[];
for k=1:size(RVELOCITYLIST,1)
        if VELINTERVAL(1)<CALVELOCITYLIST(k,3) && CALVELOCITYLIST(k,3)<VELINTERVAL(2)
        if ANGINTERVAL(1)<CALVELOCITYLIST(k,4) && CALVELOCITYLIST(k,4)<ANGINTERVAL(2)
            VPLOT=[VPLOT;CALVELOCITYLIST(k,:)];
            %[i_abs,j_abs,v,a,VELVECTOR(1),VELVECTOR(2),RMSE(i,j)];
        end
        end
end

%sort vector list
SORTV=sortrows(VPLOT,3);

%show numerical results
meanvelocity=mean(SORTV(:,3));
stdvelocity=std(SORTV(:,3));
errvelocity=stdvelocity/sqrt(size(SORTV,1));
meanangle=mean(SORTV(:,4));
stdangle=std(SORTV(:,4));
errangle=stdangle/sqrt(size(SORTV,1));
SIGNALFRAME=SIGNAL(FrameYLim(1):FrameYLim(end),FrameXLim(1):FrameXLim(end));
framearea=floor(bwarea(SIGNALFRAME))*pixelcalfactor^2;%estimate of tissue area in frame
area=floor(bwarea(SIGNAL))*pixelcalfactor^2;%estimate total excitable tissue area

fprintf(['mean velocity of selected vectors: (',num2str(meanvelocity),' +- ',num2str(errvelocity),') mm/s \n']);
fprintf(['mean angle: (',num2str(meanangle),' +- ',num2str(errangle),') degrees \n']);
fprintf(['window size (w x h) :  ',num2str(frameumwidth),' x ',num2str(frameumheight),' um \n']);
fprintf(['excitable area in window: ',num2str(framearea),' um^2 \n']);
fprintf(['total excitable area:',num2str(area),' um^2 \n']);
%% plot velocity vectors
%SYNTAX: VPLOT=[VPLOT;[i,j,i_abs,j_abs,vel,ang,Vi(i,j),Vj(i,j)];
if plotvectors==1
    figure(heartfigurezoom)
    hold on
    for k=1:size(SORTV,1)
        %SORTV=[i_abs,j_abs,v,a,VELVECTOR(1),VELVECTOR(2),RMSE(i,j)];
        x=SORTV(k,2);y=SORTV(k,1);Vx=SORTV(k,6)*vectorscale;Vy=SORTV(k,5)*vectorscale;
        %plot selected velocity vectors
        quiver(x,y,Vx,Vy,0,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth+0.25);
    end
    hold off
end
%% plot sum vector
if plotsumvector==1;%plot sum velocity vector
    figure(heartfigure);hold on
    SUMVECTOR=sum([SORTV(:,5),SORTV(:,6)])*sumvectorscale;
    if normsumvector>0
        SUMVECTOR=SUMVECTOR/norm(SUMVECTOR,2)*normsumvector;
    end
    quiver(FRAMECENTER(2)-round(SUMVECTOR(2)/2),FRAMECENTER(1)-round(SUMVECTOR(1)/2),SUMVECTOR(2),SUMVECTOR(1),0,...
        'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
    hold off
end
%save sum vector
if exist('NSUMVECTOR')==1 && sframe>0
    asksavevel=input('save sum vector(1=yes)?');
    if asksavevel==1
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
basefile=[stackfile(1:end-4),'-',description,'-delay',num2str(lowdelay),'-range',num2str(range),'-WIEN',num2str(SPATWIEN(1)),'-MED',num2str(SPATMED(1))];

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
    EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
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
    EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
    %write image to file
    p=[imagepath,basefile,'-ZOOMumxy',num2str(round(frameumwidth)),'x',num2str(round(frameumheight)),'.tif'];
    imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
    if vvectors==1
        p=[imagepath,basefile,'-VELOCITIES.mat'];
        %save(p,'VELOCITYLIST','VSTRING','ZRMSELIST','VPLOT','VRMSEPLOT','ax','ay','calvfactor_i','calvfactor_j','pixelcalfactor')
    end

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