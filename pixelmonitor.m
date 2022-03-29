P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
windowsize=0.75.*[maxwindowsize,maxwindowsize];%[width,height]*screensize
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/2);%distance of windows from bottom of screen
sep=round((screenwidth-2*windowsize(1)-2*offsetx)/2);%distance between windows
resizefactor=10;%resize divergence image for plotting
outputresolution=[1024,1024];%output image resolution [height length]
%% delete signals outside the tissue area caused by light scattering
SIGNALS=ASIGNALPIXELS;%combines ASIGNALPIXELS and SCATTER arrays
SIGNALS(SCATTER==0)=0;
DATA=NORMDATA;
for i=1:size(SIGNALS,1)
    for j=1:size(SIGNALS,2)
        if SIGNALS(i,j)==0
            DATA(i,j,:)=0;      
        end
     end
end    
%% define colors

%choose background color
%bcolor=input('background color (w) white (k) black (r) red (g) green (b) blue >','s');
bcolor='w';

%fix frame color
framecolor='r';

if bcolor=='w'
   brgb=[1 1 1];
   tx='k';
elseif bcolor=='k'
    brgb=[0 0 0];
    tx='w';
elseif bcolor=='r'
    brgb=[1 0 0];
    tx='w';
elseif bcolor=='g'
    brgb=[0 1 0];
    tx='k';
elseif bcolor=='b'
    brgb=[0 0 1];
    tx='w';
end
    


%define axes for picture display

heartfigure = figure('Name','VOLTAGE PLOT','MenuBar','none','Units','pixels','Position',[offsetx,offsety,windowsize],'Color','w','Visible','on');

%set axes for figheart
figure(heartfigure);
heart=axes('Position',[0.2 .3 .7 .7],'Visible','on','Drawmode','fast');
annotation=axes('Position',[0 .4 .4 .8],'Color','w','Visible','off','Drawmode','normal');

%set axes for figtrace
trace=axes('Position',[0.1 0.05 0.8 0.2],'Color','w','Visible','on','Drawmode','normal');


%display frames

frame=1;
button=0;
colormap jet(65535);

%choose frame
while isempty(frame)==0    
    
    %plot frame
    I=DATA(:,:,frame);
    
    [XI,map]=gray2ind(I,256);%convert 16 bit intensity image into 16 bit indexed image 
    J=ind2rgb(XI,jet(256));%convert indexed image into truecolor image using colormap jet
    
    %stain background
    for i=1:size(SIGNALS,1)
        for j=1:size(SIGNALS,2)
            if SIGNALS(i,j)==0
                J(i,j,:)=brgb;
            end
        end
    end

    %plot heart
    axes(heart),image(J),set(heart,'TickDir','out')
   
    %plot first and last activated pixels
    if exist('numpol')==1 && exist('ACTIVATIONLIST')==1 && numpol>0
        for k=first:first+numpol
            if k==1
                rectangle('Position',[(ACTIVATIONLIST(k,2)-0.5),(ACTIVATIONLIST(k,1)-0.5),1,1],'FaceColor','g');  
            else
                rectangle('Position',[ACTIVATIONLIST(k,2)-0.5,ACTIVATIONLIST(k,1)-0.5,1,1],'LineWidth',1,'EdgeColor','g');      
            end
        end
        
        for k=last-numpol+1:size(ACTIVATIONLIST,1)
            if k==size(ACTIVATIONLIST,1)
               rectangle('Position',[(ACTIVATIONLIST(k,2)-0.5),(ACTIVATIONLIST(k,1)-0.5),1,1],'FaceColor','r');   
            else
                rectangle('Position',[ACTIVATIONLIST(k,2)-0.5,ACTIVATIONLIST(k,1)-0.5,1,1],'LineWidth',1,'EdgeColor','r');  
            end
        end
    end

    %display selected pixels if available
    if exist('N')==1 && isempty(N)==0
        axes(heart)
        for i=1:size(N,1)
            rectangle('Position',[N(i,2)-0.5,N(i,1)-0.5,1,1],'LineWidth',1,'EdgeColor','w');
        end
    end
    
    %display frame if available
    if exist('framewidth')==1 && exist('frameheight')==1
        if showframe==1
            rectangle('Position',[FrameXLim(1)-0.5,FrameYLim(1)-0.5,framewidth,frameheight],'LineWidth',1,'EdgeColor','w')
        end
    end
    
    %plot time and frame number
    clear framedisplay;
    sframe=sprintf(['%0',num2str(length(num2str(size(DATA,3)))),'.0f'],frame);
    framedisplay(1)={'frame'};
    framedisplay(2)={num2str(sframe)};
    framedisplay(3)={num2str(size(DATA,3))};
    axes(annotation),cla;
    text(0,0.65,framedisplay,'color','k','FontSize',17)

    %input frame number
    oldframe=frame;
    if exist('ACTIVATIONLIST')==1
        fprintf(['First active frame: ',num2str(round(ACTIVATIONLIST(first,3))),'\n']);
        fprintf(['Last active frame: ',num2str(round(ACTIVATIONLIST(last,3))),'\n']);
    end
    frame=input(['Current frame: ',num2str(frame),'/',num2str(size(DATA,3)),' frame>']);
end

frame=oldframe;%remember last frame

%select pixels in current frame
if exist('N')==0 || isempty(N)==1
    button=0;
    figure(heartfigure);axes(heart);
    N=[];
    P=[40,40];%default pixel

    while button~=27
    
        %plot trace of selected pixel
        if ASIGNALPIXELS(P(1),P(2))==1
            PIXEL=squeeze(DATA(P(1),P(2),:));
        else
            PIXEL=-double(squeeze(BLEACHDATA(P(1),P(2),:)));
            %re-scale PIXEL
            PIXEL=PIXEL-mean(PIXEL(1:500));
            [pval,pidx]=max(PIXEL);
            PIXEL=PIXEL/mean(PIXEL(pidx-50:pidx+50));
        end
        axes(trace)
        plot(PIXEL)
    
        %plot previously selected pixels
        axes(heart)
        for i=1:size(N,1)
        rectangle('Position',[N(i,2)-0.5,N(i,1)-0.5,1,1],'LineWidth',1,'EdgeColor','w');
        end
    
        PreviousP=P;
        xi=[];yi=[];
        [xi,yi,button]=ginput(1);
    
        if isempty(button)==1 %if enter was pressed, P will be empty
            P=PreviousP;
            N=[N;P];
            button=0;
        elseif isempty(button)==0 && button~=27
            if ASIGNALPIXELS(round(yi),round(xi))==1
                P=[round(yi),round(xi)];
        else
            if exist('BLEACHDATA')==1
                P=[round(yi),round(xi)];
            else
                fprintf('Load BLEACHDATA array! \n')
            end      
        end
    end
    
    
    %if enter is pressed, save selected pixel
    if isempty(button)==1
        N=[N;previousP];
        button=0;
    end
    end

    
end
%% image
pause;
askexport=0;
askexport=input('[RETURN] export frame as tif \n [1] move SIGNALPIXELS <-> NOSIGNALPIXELS');
if isempty(askexport)==1    
%% create overlay image

    %prepare frame for image export
    mapfig = figure('Name','CONTOUR PLOT ZOOM','MenuBar','none','Units','pixels','Position',[offsetx+sep+windowsize(1),offsety,windowsize],'Color','w','Visible','on');
    mapaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);
    
    %set image resizefactor
    if exist('fresizefactor')==1
        resizefactor=fresizefactor;%fluorescence resize factor from imregist
    else
        resizefactor=10;
    end
    
    %save frame
    MAP=squeeze(NORMDATA(:,:,frame));
    %resize frame and SIGNAL image
    RMAP=imresize(MAP,resizefactor,'bicubic');
    RMAP(RMAP<0)=0;RMAP(RMAP>1)=1;
    RSIGNALS=imresize(SIGNALS,resizefactor,'bicubic');
    RSIGNALS(RSIGNALS<0.5)=0;RSIGNALS(RSIGNALS>=0.5)=1;    
    
    %convert to RGB
    [IMAP,colmap]=gray2ind(RMAP,1024);
    RGBMAP=ind2rgb(IMAP,jet(1024));
    
    %color background white
        for i=1:size(RGBMAP,1)
            for j=1:size(RGBMAP,2)
                if RSIGNALS(i,j)==0
                    RGBMAP(i,j,:)=[1,1,1];
                end
            end
        end
        
    %plot image
    axes(mapaxes)
    image(RGBMAP);
    set(mapaxes,'TickLength',[0 0])
    set(mapfig,'visible','on');
    set(mapaxes,'YDIR','reverse');
    set(mapaxes,'DataAspectRatio',[1 1 1]);

    %add recatangles
    if exist('N')==1
        for i=1:size(N,1)
            rectangle('Position',[(N(i,2)-0.5)*resizefactor,(N(i,1)-0.5)*resizefactor,resizefactor,resizefactor],'LineWidth',1.5,'EdgeColor','w');
        end
    end
    
    %plot first and last activated pixels
    if exist('ACTIVATIONLIST')==1 && numpol>0
        rectangle('Position',[(ACTIVATIONLIST(first,2)-0.5)*resizefactor,(ACTIVATIONLIST(first,1)-0.5)*resizefactor,resizefactor,resizefactor],'FaceColor','g');  
        rectangle('Position',[(ACTIVATIONLIST(last,2)-0.5)*resizefactor,(ACTIVATIONLIST(last,1)-0.5)*resizefactor,resizefactor,resizefactor],'FaceColor','r');  
    end
    
    %display frame if available
    if exist('framewidth')==1 && exist('frameheight')==1
        if showframe==1
            axes(mapaxes)
            rectangle('Position',resizefactor*[FrameXLim(1)-0.5,FrameYLim(1)-0.5,framewidth,frameheight],'LineWidth',1,'EdgeColor','w')
        end
    end
    
    %add scalebar
    scalefactor=160/72;%[um/pixels, calibration from 7/26/2007]
    unitlength=20;%length of scalebar in um
    scalebarlength=unitlength/scalefactor;%length of scalebar in pixels (80x80 image)
    fractionscalebar=0.10;
    unitstring=[num2str(unitlength),' \mum'];
    scalebarxdata=resizefactor*(80-scalebarlength+1):resizefactor*fractionscalebar:resizefactor*80;
    scalebarydata=resizefactor*80*ones(1,length(scalebarxdata));
    %draw scalebars at the bottom right corner of image
    xoffset=-40;
    yoffset=-20;
    line('xdata',scalebarxdata+xoffset,'ydata',scalebarydata+yoffset,'color','k','LineWidth',2);
    line('xdata',scalebarydata+xoffset,'ydata',scalebarxdata+yoffset,'color','k','LineWidth',2);
    %add unit length to horizontal scalebar
    text(resizefactor*70+xoffset+35,resizefactor*78+yoffset,unitstring,'color','k','FontWeight','bold')
%% export image
startpath=datapath(1:end-1);

imagepath=uigetdir(startpath,'SELECT FOLDER FOR IMAGE'); 
imagepath=[imagepath,'\'];
description=input('DESCRIPTION>','s');
basefile=['VOLTAGEMAP-FRAME',num2str(frame),'-',description];


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

%% process pixels
    
    if isempty(N)==0
        %collect pixeldata in columns and save in image folder
        PIXELDATA=[];
        for i=1:size(N,1)
        PIXELDATA(:,i)=squeeze(DATA(N(i,1),N(i,2),:));
        end
    
        %POST-PROCESSING DATA
        PIXELSMATRIX=N;%save pixel coordinates in matrix space
        %correct offset and re-normalize
        PIXELDATAPOS=[];%same data, re-aligned with maximum at maxscans/2
        idx=[];
        maxscans=1700;
        beforemaxscans=850;%850 scans = 425ms
        aftermaxscans=2000;
    
        %find areas to correct offset and interval to search for maximum
        figure,plot(PIXELDATA);
        fprintf('define offset interval \n');
        [OX,OY,B]=ginput(2);
    
        if exist('UPSTROKERANGE')==1
            searchmax=[UPSTROKERANGE(1),UPSTROKERANGE(end)];
        else
            searchmax=[1,size(PIXELDATA,1)];%search for local maximum in [delay:delay+searchmax]
        end
            
        XUPSTROKE=[];YUPSTROKE=[];SUPSTROKE=[];DSUPSTROKE=[];PIXELUPSTROKE=[];
        for i=1:size(PIXELDATA,2)
            % extract selected pixel and calculate derivative
            PIXELDATA(:,i)=PIXELDATA(:,i)-mean(PIXELDATA(round(OX(1)):round(OX(2)),i));
            %re-normalize
            delay=round(OX(2));
            [val,idx]=max(PIXELDATA([searchmax(1):searchmax(2)],i));idx=idx+searchmax(1)-1;
            PIXELDATA(:,i)=PIXELDATA(:,i)/val;
            
            %cut out upstroke and calculate derivative
            %select baseline of fit
            lthreshold=0.2;%APA to define start of upstroke
            startupstroke=idx;while PIXELDATA(startupstroke,i)>lthreshold && startupstroke>1,startupstroke=startupstroke-1;end
            
            %cut out upstroke and calculate derivative
            totalscans=500;
            moredataleft=round(0.05*scanrate);%take an additional 50 scans to the left of startupstroke
            moredataright=totalscans-idx+startupstroke-moredataleft-1;
            XUPSTROKE(:,i)=[startupstroke-moredataleft:idx+moredataright];

            YUPSTROKE(:,i)=PIXELDATA(XUPSTROKE(:,i),i);
            YUPSTROKE(:,i)=YUPSTROKE(:,i)-mean(YUPSTROKE(1:10,i));
            SUPSTROKE(:,i)=smooth(YUPSTROKE(:,i),3,'moving');
            DSUPSTROKE(:,i)=diff3(SUPSTROKE(:,i));
            %DSUPSTROKE(:,i)=DSUPSTROKE(:,i)/max(DSUPSTROKE(:,i));
            
            %save data in single array
            PIXELUPSTROKE=[PIXELUPSTROKE,XUPSTROKE(:,i),YUPSTROKE(:,i),SUPSTROKE(:,i),DSUPSTROKE(:,i)];
            
            %position data with maximum at beforemaxscans
            %find 90% of maximum
            [val,maxidx]=max(PIXELDATA(searchmax(1):searchmax(2),i));
            maxidx=maxidx+searchmax(1)-1;
            maxthreshold=0.9;
            edge=maxidx;while PIXELDATA(edge,i)>maxthreshold;edge=edge-1;end;edge=edge+1;
            if (edge<beforemaxscans && (size(PIXELDATA,1)-edge)<aftermaxscans)
                PIXELDATAPOS(:,i)=[zeros(beforemaxscans-edge,1);PIXELDATA(:,i);zeros(aftermaxscans-(size(PIXELDATA,1)-edge),1)];
    
            elseif (edge<beforemaxscans && (aftermaxscans<(size(PIXELDATA,1)-edge)))
                PIXELDATAPOS(:,i)=[zeros(beforemaxscans-edge,1);PIXELDATA(1:edge+aftermaxscans,i)];
    
            elseif (beforemaxscans<edge && (size(PIXELDATA,1)-edge)<aftermaxscans)
                PIXELDATAPOS(:,i)=[PIXELDATA(edge-beforemaxscans+1:end,i);zeros(aftermaxscans-(size(PIXELDATA,1)-edge),1)];
    
            elseif (beforemaxscans<=edge && aftermaxscans<=(size(PIXELDATA,1)-edge))
                PIXELDATAPOS(:,i)=PIXELDATA(edge-beforemaxscans+1:edge+aftermaxscans,i);
            end
            
        end
        


% save pixel data
        
    pixelmat=[imagefile(1:end-4),'.mat'];
    pixeltxtdata=[imagefile(1:end-4),'-time.txt'];
    pixeltxtmax=[imagefile(1:end-4),'-max.txt'];
    pixeltxtmaxmean=[imagefile(1:end-4),'-mean.txt'];
    pixelup=[imagefile(1:end-4),'-upstroke.txt'];
    pixelupmean=[imagefile(1:end-4),'-meanup.txt'];
    pixeldiffup=[imagefile(1:end-4),'-dup.txt'];
    pixeldiffup=[imagefile(1:end-4),'-meandup.txt'];
      
    save([imagepath,pixelmat],'PIXELDATA','PIXELDATAPOS','stackfile','PIXELSMATRIX','XUPSTROKE',...
        'YUPSTROKE','SUPSTROKE','DSUPSTROKE');
    save([imagepath,pixeltxtdata],'PIXELDATA','-ascii');
    save([imagepath,pixeltxtmax],'PIXELDATAPOS','-ascii');
    save([imagepath,pixelup],'SUPSTROKE','-ascii');
    save([imagepath,pixeldiffup],'DSUPSTROKE','-ascii');
    

% 
    end
%%
    end
% transfer pixels between SIGNALPIXELS and NOSIGNALPIXELS if askexport==1  
    %transfer chosen pixels
    if isempty(askexport)==0 && askexport==1
%% transfer pixels
        for m=1:size(N,1)
            i=N(m,1);j=N(m,2);
            if ASIGNALPIXELS(i,j)==0
                ASIGNALPIXELS(i,j)=1;
            else 
                ASIGNALPIXELS(i,j)=0;
            end
        end
                
%% end of if statement
    end
    