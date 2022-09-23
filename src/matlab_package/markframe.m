% FRAMECENTER=DATABASE(1).v_windowcenter;
% setframeumwidth=round(DATABASE(2).v_windowx);
% setframeumheight=round(DATABASE(2).v_windowy);
showframe=1;
mapdata=2;%0=CA MAP 1=APD MAP 2=ISOCHRONAL MAP
FRAMELOCATIONS=NFRAME;
for i=1:size(FRAMELOCATIONS,1)-1
FrameXLim=FRAMELOCATIONS(i,1:2);
FrameYLim=FRAMELOCATIONS(i,3:4);
framewidth=FrameXLim(2)-FrameXLim(1);
frameheight=FrameYLim(2)-FrameYLim(1);
FRAMECENTER=[mean(FrameYLim),mean(FrameXLim)];

%calculate correct frameumwidth
frameumwidth=framewidth*pixelcalfactor;frameumheight=frameheight*pixelcalfactor;

% add frame to mark selected area

if showframe==1
    if exist('RLOW')==1 && mapdata==0
        axes(mapaxes);%CA MAP
        rectangle('Position',resizefactor*[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
    
    elseif exist('APDMATRIX')==1 && mapdata==1
        axes(mapaxes);%APD MAP
        rectangle('Position',resizefactor*[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--')
    
    else
        axes(isoc);%voltage data
        rectangle('Position',[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--');
    end
    
    
end
    
if exist('RGBSHIFTB')==1
        axes(heartimageaxes)
        rectangle('Position',fresizefactor*[FrameXLim(1),FrameYLim(1),framewidth,frameheight],'LineWidth',frameweight,'EdgeColor',framecolor,'LineStyle','--');
end

%plot sum vector
%normsumvector=10;%normalizes sum vector so that all sumvectors have the same length
if exist('NSUMVECTOR')==1 && mapdata==2;%plot sum velocity vector
    axes(isoc);hold on
    SUMVECTOR=NSUMVECTOR(i,:);
    if normsumvector>0
        SUMVECTOR=SUMVECTOR/norm(SUMVECTOR,2)*normsumvector;
    end
    quiver(FRAMECENTER(2)-round(SUMVECTOR(2)/2),FRAMECENTER(1)-round(SUMVECTOR(1)/2),SUMVECTOR(2),SUMVECTOR(1),0,...
        'LineWidth',sumvectorwidth,'Color',sumvectorcolor,'MaxHeadSize',1);
    hold off
end
pause
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

%mapdata=0:CA mapdata=1:APD mapdata=2:ISOCHRONAL MAP
if mapdata==0 || mapdata==1
    basefile=[stackfile(1:end-4),'-map-',description];
else
    basefile=[stackfile(1:end-4),'-',description,'-delay',num2str(lowdelay),'-range',num2str(range),'-WIEN',num2str(SPATWIEN(1)),'-MED',num2str(SPATMED(1))];
end

if isempty(askexport)==1
    
    %WRITE CONTOUR PLOT
    if mapdata==0 || mapdata==1
        figure(mapfig)
    else   
        figure(heartfigure)
    end
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
    if mapdata==2
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
    
    %WRITE HEART IMAGE
    if exist('RGBSHIFTB')==1
        axes(heartimageaxes)
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
end