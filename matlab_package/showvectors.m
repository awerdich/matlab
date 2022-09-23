%plot velocity vectors
%showcontours has to be run before
%% parameters
showframe=1;
showcontourlines=0;%0:do not show contour lines
fillcontour='off';
lineweight=2.5;%line width of contour plot
lineweightzoom=0.3;%line width of zoom contour plot
vectorscale=0.3;%1.5 for 24hpf
normvectorscale=0.5;%plot normalized vectors 0: do not normalize [0.5]
maxheadsize=4;
vectorwidth=1.5;%[pts] vector line width [2.5]
vectorcolor=[0 0 0];
%TIME LIMITS
lowdelay=-9;%time delay to plot after smallest value
range=500;
dt=10;
lowlimit=STIMES(1,1)+lowdelay;
highlimit=lowlimit+range;
COLORRANGE=[lowlimit,highlimit];
division=(highlimit-lowlimit)/2;
tickdivision=numcolors/(highlimit-lowlimit)*division;

%define window size based on actual screen size
P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/5);%distance of windows from bottom of screen
sep=round((screenwidth-2*windowsize(1)-2*offsetx)/2);%distance between windows
%aspect ratios and output resolution
matrixaspectratio=size(C,2)/size(C,1);%aspectratio of data input matrix [width/height]
outputresolution=[1024,round(1024/matrixaspectratio)];%output image resolution [width, height]
frameaspectratio=framewidth/frameheight;
if maxwindowsize<windowsize(2)/frameaspectratio
    wscalefactor=0.75*windowsize(2)/maxwindowsize/frameaspectratio;
else
    wscalefactor=0.75;
end
windowsize=wscalefactor*[maxwindowsize,maxwindowsize];%[width,height]*screensize   
    
outputresolutionzoom=[outputresolution(1),round(outputresolution(2)/frameaspectratio)];%resolution of zoom image [width,height]
%% prepare figures and axes

%define axes for picture display
heartfigure = figure('Name','CONTOUR PLOT','MenuBar','none','Units','pixels','Position',[offsetx,offsety,windowsize(1),windowsize(2)/matrixaspectratio],'Color','w','Visible','on');
heartfigurezoom = figure('Name','CONTOUR PLOT ZOOM','MenuBar','none','Units','pixels','Position',[offsetx+sep+windowsize(1),offsety,windowsize(1),windowsize(2)/frameaspectratio],'Color','w','Visible','on');

%set axes for figheart
figure(heartfigure);
isoc=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);

figure(heartfigurezoom);
isoczoom=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);
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
% add frame to mark selected area
axes(isoc);
if showframe==1
    set(isoc,'DrawMode','normal');
    rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
    %rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--');
end
%% define zoom contour plot
caxis(COLORRANGE);
[CONTZOOM,hzoom]=contour(isoczoom,INTERPX,INTERPY,CONTM,[lowlimit:dt:highlimit],'LineColor','k');
set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
%set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
set(hzoom,'Fill',fillcontour,'LineWidth',lineweightzoom);
%% draw polygon
if exist('POLYROI','var')==1 && isempty(POLYROI)==0
    axes(isoczoom);
    polygonhandle=impoly(isoczoom,POLYGONPOSITIONS);
end
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
VELINTERVAL=[meanvelall-4*stdvelall,meanvelall+4*stdvelall];
VPLOTBIN=zeros(size(SIGNALFRAMEINTERP_0));%binary matrix of filtered pixels
VPLOTLIST=[];%list of filtered velocities
for i=1:size(SIGNALFRAMEINTERP_0,1)
    for j=1:size(SIGNALFRAMEINTERP_0,2)
        if SIGNALFRAMEINTERP_0(i,j)==1
            if (0<RMSEFRAME(i,j) && RMSEFRAME(i,j)<maxrmse)
                v=norm([CALVFRAMEi(i,j);CALVFRAMEj(i,j)],2);
                if VELINTERVAL(1)<v && v<VELINTERVAL(2)
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
%SYNTAX: VPLOT=%[i_abs,j_abs,v,a,VELVECTOR(1),VELVECTOR(2),RMSE(i,j)];
%plotvectors=0;
if plotvectors==1
    VECTORLIST=[];
    for i=1:size(VPLOTBIN,1)
        for j=1:size(VPLOTBIN,2)
            if VPLOTBIN(i,j)==1
                x=XFRAMEINTERP_0(1,j);
                y=YFRAMEINTERP_0(i,1);
                vx=CALVFRAMEj(i,j);
                vy=CALVFRAMEi(i,j);
                if normvectorscale>0
                    %normalize velocities
                    v=norm([vx,vy],2);
                    vx=vx/v*normvectorscale;
                    vy=vy/v*normvectorscale;
                    vectorscale=0;%set vectorscale to 0 to use normvectorscale in plot function
                end
                if exist('POLYROI','var')==1 && isempty(POLYROI)==0
                    if POLYROI(i,j)==1
                    %plot vector [vx,vy] at position [x,y]
                    VECTORLIST=[VECTORLIST;[x,y,vx,vy]];
                    %quiver(x,y,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
                    end
                else
                    %plot vector [vx,vy] at position [x,y]
                    VECTORLIST=[VECTORLIST;[x,y,vx,vy]];
                    %quiver(x,y,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
                end
            end
        end
    end
% plot vectorsfigure(heartfigurezoom)
axes(isoczoom)
hold on
for i=1:1:length(VECTORLIST)
    x=VECTORLIST(i,1);y=VECTORLIST(i,2);vx=VECTORLIST(i,3);vy=VECTORLIST(i,4);
    quiver(x,y,vx,vy,vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
end
hold off
end
%remove contour lines, if requested
if showcontourlines==0
    set(hzoom,'Visible','off');
end
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
    EXPORTIMAGE=imresize(IMAGEC,[outputresolutionzoom(2),outputresolutionzoom(1)],'bicubic');
    %write image to file
    p=[imagepath,basefile,'-ZOOMumxy',num2str(round(frameumwidth)),'x',num2str(round(frameumheight)),'.tif'];
    imwrite(EXPORTIMAGE,p,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
    if vvectors==1
        p=[imagepath,basefile,'-VELOCITIES.mat'];
        save(p,'scanrate','pixelcalfactor_x','sumvectorscale','vectorscale','normsumvector','sumvectorwidth','vectorwidth',...
            'vectorcolor','sumvectorcolor','VELINTERVAL','VPLOTBIN','FRAMECENTER');
    end
end