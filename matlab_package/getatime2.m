%get activation times from database
%% field definitions
region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');

if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
|| strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
|| strcmp(region,'vic')==1 || strcmp(region,'voc')==1

regionpol=[region,'_POL'];

end
%% LOOP
% find activation times for pixels that are far away
WAVEFRONTS=[];SET=[];
for id=1:length(DATABASE)
    %% extract activation map from database
    scanrate=DATABASE(id).scanrate;%scans/s
    if isfield(DATABASE,'pixelcalfactorj')==1
        pixelcalfactor=DATABASE(id).pixelcalfactorj;%[um/pixel]
    elseif isfield(DATABASE,'pixecalfactor')==1
        pixelcalfactor=DATABASE(id).pixelcalfactor;
    end
    if isempty(DATABASE(id).(regionpol))==0
        SET=[SET;id];%dataset used
    %activation matrix
    POL=DATABASE(id).(regionpol);
    %% normalize 
    %determine data range
    T=[];
    for i=1:size(POL,1)
        for j=1:size(POL,2)
            if POL(i,j)>0
                T=[T;POL(i,j)];
            end
        end
    end
    tmin=min(T);tmax=max(T);%activation time data limit
    lowlimit=0;highlimit=1.0;%additional limit for lowest and highest activation times
    
    %re-map activation times
    NPOL=zeros(size(POL));
    LNPOL=zeros(size(POL));
    for i=1:size(POL,1)
        for j=1:size(POL,2)
            if POL(i,j)>0
                    %map entire activation time into [0,1]
                    NPOL(i,j)=(POL(i,j)-tmin)/(tmax-tmin);
                    %map activation times into [lowlimit,highlimit]
                    if NPOL(i,j)<lowlimit
                        LNPOL(i,j)=0;
                    elseif lowlimit<=NPOL(i,j) && NPOL(i,j)<highlimit
                        LNPOL(i,j)=(NPOL(i,j)-lowlimit)/(highlimit-lowlimit);
                    else
                        LNPOL(i,j)=1;
                    end
            end
        end
    end
    
    %Re-size data matrix
    resizefactor=2;%resize divergence image for plotting
    %signal matrix
    SIGNAL=zeros(size(POL));
    SIGNAL(POL>0)=1;
    %re-size data
    RNPOL=imresize(LNPOL,resizefactor,'bicubic');
    RSIGNAL=imresize(SIGNAL,resizefactor,'bicubic');
    RSIGNAL(RSIGNAL<0.9)=0;RSIGNAL(RSIGNAL>=0.9)=1;
    %fix interpolated activation matrix at signalk boundary
    RNPOL(RSIGNAL==0)=0;
    
    %thresholding 
    dt=0.1;%time interval for thresholding
    t=dt;%start time
    n=floor(1/dt);%number of time steps
    TPOL=zeros(size(RNPOL,1),size(RNPOL,2),n);
    ISLANDS=[];%number of wavefronts
    while t+dt<=1
        TINTERVAL=[t;t+dt];
        TPOL=zeros(size(RNPOL));
        for i=1:size(RNPOL,1)
            for j=1:size(RNPOL,2)
                %exclude pixels that are off
                if RSIGNAL(i,j)==1
                    if TINTERVAL(1)<=RNPOL(i,j) && RNPOL(i,j)<=TINTERVAL(2)
                        TPOL(i,j)=1;
                    end
                end
            end
        end
        %calculate the connectivity of TPOL
        CC=bwconncomp(TPOL,8);
        ISLANDS=[ISLANDS;[id,t,CC.NumObjects]];%islands as 
        t=t+dt;
    end
    %maximum number of islands
    [val,idx]=max(ISLANDS(:,3));
    WAVEFRONTS=[WAVEFRONTS;ISLANDS(idx,:)];
    save WAVEFRONTS.txt WAVEFRONTS -ascii -tabs
    end
    end
%% plot data    
P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
windowsize=0.75.*[maxwindowsize,maxwindowsize];%[width,height]*screensize
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/2);%distance of windows from bottom of screen
sep=round((screenwidth-2*windowsize(1)-2*offsetx)/2);%distance between windows
outputresolution=[1024,1024];%output image resolution [height length]    
mapfig = figure('Name','ACTIVATION MAP','MenuBar','none','Units','pixels','Position',[offsetx+sep+windowsize(1),offsety,windowsize],'Color','w','Visible','on');
mapaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','normal','YDir','reverse','TickLength',[0 0]);

%convert image to RGB
[IMAP,colmap]=gray2ind(mat2gray(RNPOL),1024);
RGB_RNPOL=ind2rgb(IMAP,jet(1024));
 
%change color of unmapped tissue white
for i=1:size(RGB_RNPOL,1)
    for j=1:size(RGB_RNPOL,2)
        if RSIGNAL(i,j)==0
            RGB_RNPOL(i,j,:)=[1,1,1];
        end
    end
end

%plot image
image(RGB_RNPOL)

%mark activated pixels at time t
hold on;
for i=1:size(TPOL,1)
    for j=1:size(TPOL,2)
        if TPOL(i,j)==1
            rectangle('Position',[j-0.5,i-0.5,1,1],'EdgeColor',[0,0,0],'LineWidth',3)
        end
    end
end

%save image
if exist('imagepath','var')==0
    imagepath=[];
end
imagepath=uigetdir(imagepath,'pick path');imagepath=[imagepath,'\'];
filemap=[DATABASE(id).name,'-t',num2str(t*100),'-',region,'.tif'];
imagefile=[imagepath,filemap];
axes(mapaxes)
set(mapaxes,'DataAspectratioMode','auto','TickLength',[0,0]);
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%crop image to remove black border
IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,2)-4,size(IMAGE,1)-4]);
%IMAGEC=IMAGE;
EXPORTIMAGE=imresize(IMAGEC,outputresolution,'bicubic');
%write image to file
imwrite(EXPORTIMAGE,imagefile,'tif','Compression','none','ColorSpace','rgb','Resolution',600);
set(mapaxes,'DataASpectratioMode','manual');

