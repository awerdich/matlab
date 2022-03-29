%Draw whole heart vectormap
%Assumes showcontours was run
vectorcolor=[0 0 0];
vectorwidth=1.25;%[pts] vector line width
vectorscale=0.2;%scaling factor of velocity vector field (usually 0.2 72hpf, 3 24 hpf)
lineweightzoom=0.1;%line width of zoom contour plot
frameaspectratio=(FrameXLim(2)-FrameXLim(1))/(FrameYLim(2)-FrameYLim(1));
outputresolution=[1024,1024*frameaspectratio];%output image resolution [height length]
%% filter vectors
%VELOCITYLIST=[YFRAMEINTERP_0(i,1),XFRAMEINTERP_0(1,j),v,alpha,VELVECTOR(2),VELVECTOR(1),RMSEFRAME(i,j),i,j]];%matrix coordinates 
VELINTERVAL=[meanvelall-4*stdvelall,meanvelall+2*stdvelall];
ANGINTERVAL=[-180,180];
%ANGINTERVAL=[meanangleall-stdangleall,meanangleall+stdangleall];
%ANGINTERVAL=[meanangleall-2*stdangleall,meanangleall+2*stdangleall];
%ANGINTERVAL=[meanangleall-3*stdangleall,meanangleall+3*stdangleall];
VPLOT=[];%list of filtered velocities
VPLOTBIN=zeros(size(RMSEFRAME));%binary FRAME image of pixels used to calulate velocities
for k=1:size(CALVELOCITYLIST,1)
        i=CALVELOCITYLIST(k,8);j=CALVELOCITYLIST(k,9);
        if VELINTERVAL(1)<CALVELOCITYLIST(k,3) && CALVELOCITYLIST(k,3)<VELINTERVAL(2)
        if ANGINTERVAL(1)<CALVELOCITYLIST(k,4) && CALVELOCITYLIST(k,4)<ANGINTERVAL(2)
            VPLOT=[VPLOT;CALVELOCITYLIST(k,:)];
            %[YFRAMEINTERP_0(i,1),XFRAMEINTERP_0(1,j),v,alpha,VELVECTOR(2),VELVECTOR(1),RMSEFRAME(i,j),i,j]]
            
            VPLOTBIN(i,j)=1;
        end
        end
end

%sort vector list
SORTV=sortrows(VPLOT,3);
%% define zoom contour plot
caxis(COLORRANGE);
[CONTZOOM,hzoom]=contour(isoczoom,INTERPX,INTERPY,CONTM,[lowlimit:dt:highlimit],'LineColor','k');
set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
%set(isoczoom,'YDir','reverse','Color','none','TickDir','out','XLim',FrameXLim,'YLim',FrameYLim);
set(hzoom,'Fill','off','LineWidth',lineweightzoom);
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
    % ask to change minlength
    askminlength=input(['minimum contour length [',num2str(minlength),']>']);
    %askminlength=[];
end
%% plot velocity vectors
%SYNTAX: VPLOT=%[i_abs,j_abs,v,a,VELVECTOR(1),VELVECTOR(2),RMSE(i,j)];
%plotvectors=0;
if plotvectors==1
    figure(heartfigurezoom)
    hold on
    for k=1:2:size(SORTV,1)
        %SORTV=[i_abs,j_abs,v,a,VELVECTOR(i),VELVECTOR(j),RMSE(i,j)];
        x=SORTV(k,2);y=SORTV(k,1);Vx=SORTV(k,6)*vectorscale;Vy=SORTV(k,5)*vectorscale;
        %plot selected velocity vectors
        %quiver(x,y,Vx,Vy,0,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth+0.25);
        quiver(x,y,Vx,Vy,0,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth+0.5);
    end
    hold off
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
        save(p,'scanrate','pixelcalfactor','sumvectorscale','vectorscale','normsumvector','sumvectorwidth','vectorwidth',...
            'vectorcolor','sumvectorcolor','VELINTERVAL','ANGINTERVAL','VPLOT','SORTV','SUMVECTOR','FRAMECENTER');
    end
end