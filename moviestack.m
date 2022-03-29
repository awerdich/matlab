[DATA,SIGNAL,SF,scanrate,NSTIM]=moviestackprep(NORMDATA,ASIGNALPIXELS,SCATTER,STIM);
opengl('software');
resizefactor=10;
firstframe=1;
lastframe=size(DATA,3);
newpixel=1;
traceylimit=[-0.2,1.4];%limit for Y axis on trace
lright=[0.75 0.22 0.05 0.2];%right side colorbar 'position'
lleft=[0.01 0.22 0.05 0.2];%left side colorbar 'position'
uleft=[0.01 0.75 0.05 0.2];
uright=[0.9 0.75 0.05 0.2];
colorbarposition=lleft;
if exist('ACTIVATIONLIST')==1
    first=10;last=size(ACTIVATIONLIST,1);
end

overlaythreshold=0.2;%fluorescence display threshold
smallmovie=1;%0=only large movie, 1=small and large movie

%PIXEL CALIBRATION FACTOR
pixelcalfactor=pixelcalfactor_x;
unitlength=20;%length of scalebar in um
scalebarlength=unitlength/pixelcalfactor;%length of scalebar in pixels (80x80 image)

%% create figure and axes
aspectwh=size(DATA,2)/size(DATA,1);

%figure and axes for large frame
fig = figure('Name','Propagating Action Potentials','MenuBar','none','Units','pixels','Position',[20 100 640 640/aspectwh],'Color','k','Visible','on');
set(gcf,'DoubleBuffer','on'); 
heart=axes('Position',[0 .19 1 .82],'Visible','off','Drawmode','fast');
set(heart,'TickDir','out');
trace=axes('Position',[0 0.01 1 .17],'Color','w','Visible','on','Drawmode','fast','YLim',traceylimit);

%figure and axes for small frame
if smallmovie>0
    fig2= figure('Name','Propagating Action Potentials','MenuBar','none','Units','pixels','Position',[680 380 320 320/aspectwh],'Color','k','Visible','on');
    heart2=axes('Position',[0 0 1 1],'Visible','off','Drawmode','fast');
    set(heart2,'TickDir','out','TickLength',[0 0]);
end

% display frame
FIRSTFRAME=imadjust(mat2gray(squeeze(DATA(:,:,1))));
axes(heart)
imagesc(FIRSTFRAME);colormap(gray)
%% pick pixel to display
figure(fig);
axes(heart);

if exist('pixelx')==1 && newpixel==0
    button=[];
else
    button=1;
end
    
while isempty(button)==0
[xi,yi,button]=ginput(1);
if isempty(button)==0 %save pixel values if enter is pressed
pixelx=round(xi);pixely=round(yi);
end
axes(trace),plot(squeeze(DATA(pixely,pixelx,:)),'w','LineWidth',2);
set(trace,'Color','k','YLim',traceylimit)
end

%% get movie filename and path

%pick filename and movieparameter
movpath=uigetdir(datapath,'Pick folder to save movie file.');
movpath=[movpath,'\'];
movfile=[stackfile(1:end-4),'-',num2str(firstframe),'-',num2str(lastframe),'.avi'];
movfile2=[stackfile(1:end-4),'-S-',num2str(firstframe),'-',num2str(lastframe),'.avi'];

mov=VideoWriter([movpath,movfile],'Motion JPEG AVI');
mov2=VideoWriter([movpath,movfile2],'Motion JPEG AVI');
mov.Quality=100;mov2.Quality=100;
mov.FrameRate=30;mov2.FrameRate=30;
open(mov);open(mov2);

%% define scalebar
fractionscalebar=0.10;
unitstring=[num2str(unitlength),' \mum'];
if exist('fresizefactor')==1
    resizefactor=fresizefactor;%fluorescence resize factor from imregist
else
    resizefactor=10;%image resize before avi generation
end
scalebarxdata=resizefactor*(80-scalebarlength+1):resizefactor*fractionscalebar:resizefactor*80;
scalebarydata=resizefactor*80*ones(1,length(scalebarxdata));
%% main loop
%resize ASIGNALPIXELS
RSIGNAL=imresize(SIGNAL,resizefactor);
RSIGNAL(RSIGNAL<0.5)=0;RSIGNAL(RSIGNAL>=0.5)=1;

%resize and convert background image
if isempty(SF)==0
    %resize background image
    RSF=imadjust(mat2gray(imresize(SF,resizefactor,'bicubic')));
    [IDX,map]=gray2ind(RSF,256);
    RGBSF=ind2rgb(IDX,gray(256));
end

PIXEL=squeeze(DATA(pixely,pixelx,firstframe:lastframe));
PIXELSTIM=NSTIM(firstframe:lastframe);

halftracelength=length(PIXEL);
nframes=lastframe-firstframe+1;
% colorbar parameters
numcolor=100;
colormap(jet(numcolor));
CLim=[1,numcolor];%color limits
YLim=[overlaythreshold*numcolor,numcolor];
dc=0.2;%distance of color tickmarks
YTick=[overlaythreshold*numcolor:dc*numcolor:numcolor];
YTickLabel=[];
for k=1:length(YTick)
    YTickLabel=[YTickLabel,{num2str(YTick(k)/numcolor)}];
end

for frame=1:2:nframes
%% prepare single frame     
    %resize frame
    I=DATA(:,:,firstframe+frame-1);%get frame
    %re-set pixels of I below plotting threshold
    I(I<overlaythreshold)=0;
    RESI=imresize(I,resizefactor,'bicubic');
    %recolor fluorescence image border using NBORDER
    if exist('NBORDER')==1
        RESI(NBORDER==0)=0;
    end
    [RESX,mapRESX]=gray2ind(RESI,65536);
    RGBRESI=ind2rgb(RESX,jet(65536));
    
    %plot frame
    if isempty(SF)==0
        %single background image
        axes(heart)
        image(RGBSF);hold on 
            if smallmovie>0
                axes(heart2)
                image(RGBSF);hold on
            end
        %ALPHAMATRIX for fluorescence image
        transparency=0.7;
        ALPHAI=ones(size(RESI))*transparency;
        for i=1:size(RESI,1)
            for j=1:size(RESI,2)
                if RESI(i,j)<overlaythreshold+0.1
                    ALPHAI(i,j)=0;
                end
             end
        end
        axes(heart)
        %THIS COMMAND REVERSES YAXIS ON SOME SYSTEMS
        image(RGBRESI,'AlphaData',ALPHAI);
        %set(heart,'YDir','normal');
        hold off
        if smallmovie>0
            axes(heart2)
            image(RGBRESI,'AlphaData',ALPHAI);
            hold off
        end
    else
        %no background imagage
        for a=1:size(RGBRESI,1)
            for b=1:size(RGBRESI,2)
                if RSIGNAL(a,b)<1
                    RGBRESI(a,b,:)=[0,0,0];
                end
            end
        end
        axes(heart),image(RGBRESI);
        if smallmovie>0
            axes(heart2),image(RGBRESI);
        end
    end

    axes(heart);
    %Plot rectangle around selected pixel
    rectangle('Position',[(pixelx-0.5)*resizefactor,(pixely-0.5)*resizefactor,resizefactor,resizefactor],'LineWidth',1.5,'EdgeColor','w')
    
    %plot first and last activated pixels
%     if exist('ACTIVATIONLIST')==1 && exist('first')==1
%         axes(heart)
%         rectangle('Position',resizefactor*[ACTIVATIONLIST(first,2)-0.5,ACTIVATIONLIST(first,1)-0.5,1,1],'EdgeColor','g','FaceColor','g');  
%         rectangle('Position',resizefactor*[ACTIVATIONLIST(last,2)-0.5,ACTIVATIONLIST(last,1)-0.5,1,1],'EdgeColor','r','FaceColor','r');  
%         axes(heart2)
%         rectangle('Position',resizefactor*[ACTIVATIONLIST(first,2)-0.5,ACTIVATIONLIST(first,1)-0.5,1,1],'EdgeColor','g','FaceColor','g');  
%         rectangle('Position',resizefactor*[ACTIVATIONLIST(last,2)-0.5,ACTIVATIONLIST(last,1)-0.5,1,1],'EdgeColor','r','FaceColor','r');  
%     end
    
    %draw scalebars at the bottom right corner of image
    xoffset=-10;
    yoffset=-20;

%     axes(heart2)
%     line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','w','LineWidth',2);
%     line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','w','LineWidth',2);
    
%     axes(heart)
%     line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','w','LineWidth',2);
%     line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','w','LineWidth',2);
    
    %add unit length to horizontal scalebar
%     text(resizefactor*70+xoffset+20,resizefactor*78+yoffset,unitstring,'color','w','FontWeight','bold','FontSize',10)
    
    %prepare TRACE DATA
    %fluorescence
    TRACE=[PIXEL(1)+zeros(halftracelength-frame+1,1);PIXEL(1:end-1);PIXEL(end)+zeros(frame,1)];
    axes(trace),plot(TRACE,'w','LineWidth',2)
    hold on
    line('xdata',halftracelength*ones(10,1),'ydata',[0:0.1:0.9],'color','r','LineWidth',2);
    %plot(halftracelength,[0:.001:1.0],'.','LineWidth',2,'Color','r')
    axes(trace),plot(TRACE(1:halftracelength-frame+1),'k','LineWidth',2);
    xdata=[2*halftracelength-frame:length(TRACE)];
    plot(xdata,TRACE(halftracelength+(halftracelength-frame):end),'k','LineWidth',2);
    set(trace,'Color','k','YLim',traceylimit) 
    
    %display stimulus in same graph if exist
    if exist('PIXELSTIM')==1
        hold on
        TRACESTIM=[PIXELSTIM(1)+zeros(halftracelength-frame+1,1);PIXELSTIM(1:end-1);PIXELSTIM(end)+zeros(frame,1)];
        axes(trace),plot(TRACESTIM,'y','LineWidth',1)
        axes(trace),plot(TRACESTIM(1:halftracelength-frame+1),'k','LineWidth',1);
        plot(xdata,TRACESTIM(halftracelength+(halftracelength-frame):end),'k','LineWidth',1);
    end
    stime=sprintf('%2.1f',(frame-1)/scanrate*1e3);
    timedisplay(1) = {['Time ',num2str(stime),' ms']};
    text(length(TRACE)*(0.5-0.05),1.3,timedisplay,'color','w','FontSize',10,'FontWeight','bold');hold off
    %plot colorbar in heart axes
    axes(heart);
    colormap(jet(numcolor));
    %colorbar('EastOutside','position',colorbarposition,'Clim',CLim,'YLim',YLim,'YTick',YTick,'Box','off','YColor','w','XColor','k','YAxisLocation','right','YTickLabel',YTickLabel,'DrawMode','fast','TickLength',[0 0],'FontSize',10,'FontWeight','bold');
%% save frame in file  
    %save frame in file  
    drawnow;
    set(heart,'TickLength',[0 0],'visible','off');
    set(heart2,'TickLength',[0 0],'visible','off');
    set(trace,'TickLength',[0 0],'visible','off');
    drawnow;
    writeVideo(mov,getframe(fig));
    %mov = addframe(mov,currentframe);
    if smallmovie>0
        writeVideo(mov2,getframe(fig2));
        %mov2 = addframe(mov2,currentframe2);
    end
end
close(mov);close(mov2);    