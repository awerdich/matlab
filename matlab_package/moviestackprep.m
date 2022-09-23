function [DATA,SIGNAL,SF,scanrate,NSTIM]=moviestackprep(NORMDATA,ASIGNALPIXELS,SCATTER,STIM)
%% run function in workspace
camscanrate=2000;%camera scanrate
FrameXLim=[1,80];
FrameYLim=[1,80];

%build signal matrix
S=ASIGNALPIXELS;
S(SCATTER==0)=0;
SIGNAL=S(FrameYLim(1):FrameYLim(2),FrameXLim(1):FrameXLim(2));

firstframe=1;
lastframe=2500;

%total number of frames
movieframes=600;%10 seconds at 30 frames/s    

%crop data for ROI
if exist('FrameXLim','var')==1 && isempty(FrameXLim)==0
    %build new image stack 
    CDATA=NORMDATA(FrameYLim(1):FrameYLim(2),FrameXLim(1):FrameXLim(2),firstframe:lastframe);
else
    CDATA=NORMDATA(:,:,firstframe:lastframe);
end

%interpolate CDATA to obtain the desired number of frames (movieframes)
T=([1:size(CDATA,3)]'-1)/camscanrate*1.0e3;%original time axis
idt=(T(end)-T(1))/(movieframes-1);
IT=[T(1):idt:T(end)]';%new time axis
scanrate=(1/idt)*1.0e3;
NDATA=zeros(size(CDATA,1),size(CDATA,2),length(IT));
for i=1:size(CDATA,1)
    for j=1:size(CDATA,2)
        if SIGNAL(i,j)==1
            PIXEL=squeeze(CDATA(i,j,:));
            %re-normalize each trace
            OPIXEL=PIXEL-mean(PIXEL(1:10));
            NPIXEL=OPIXEL/max(OPIXEL);
            IPIXEL=interp1(T,NPIXEL,IT,'spline');
            NDATA(i,j,:)=IPIXEL;
        end
    end
end

%interpolate stimulus
S=STIM(firstframe:lastframe);
%re-normalize each trace
OS=S-mean(S(1:10));
NS=OS/max(OS);
IS=interp1(T,NS,IT,'spline');
NSTIM=IS;

%load background image
datapath=[];
[backgroundfile backgroundpath]=uigetfile([datapath,'*.tif'],'load background image taken with high-speed camera or cancel');    
if length(backgroundfile)>1
    SF=mat2gray(imread([backgroundpath,backgroundfile]));
    SF=imadjust(SF);
    SF=SF(FrameYLim(1):FrameYLim(2),FrameXLim(1):FrameXLim(2));   
else
    pause
    SF=[];
end
%% determine angle and flip reqiorements
%determine angles and flip requirements
figure;colormap(gray)
%prepare binary image to show orientation
F=SIGNAL;
angle=0;
previousangle=angle;
flip=0;
previousflip=flip;
FR=F;
newangle=0;
newflip=0;
while isempty(newangle)==0 && isempty(newflip)==0
    imagesc(FR)
    newangle=input(['Angle [ENTER]=',num2str(angle),'>']);
    if isempty(newangle)==0 && -180<=newangle && newangle<=180
        if newangle<0
            angle=360+newangle;
        else
            angle=newangle;
        end
        %ROTATE
        FR=imrotate(F,angle,'crop');
        imagesc(FR);
    end
    newflip=input(['Flip 1=VERTIVAL 2=HORIZONTAL 3=BOTH [ENTER]=',num2str(flip),'>']);
    if isempty(newflip)==0 && 0<=newflip && newflip<=3
        flip=newflip;
        %FLIP
        if flip==1 || flip==3
            FR=flipud(FR);
            imagesc(FR);
        end
        if flip==2 || flip==3
            FR=fliplr(FR);
            imagesc(FR);
        end
    end
end
%% rotate and flip data
% process data
RDATA=zeros(size(NDATA));
if angle>0 || flip>0
    hdl = waitbar(0,['ROTATING AND FLIPPING IMAGE STACK']);
    for frame=1:size(NDATA,3)
        %ROTATE
        ROTFRAME=imrotate(squeeze(NDATA(:,:,frame)),angle,'crop');
        ROTSF=imrotate(SF,angle,'crop');
        
        %FLIP
        if flip==1 || flip==3
            FLIPFRAME=flipud(ROTFRAME);
            FLIPSF=flipud(ROTSF);
        else
            FLIPFRAME=ROTFRAME;
            FLIPSF=ROTSF;
        end
        
        if flip==2 || flip==3
            FLIPFRAME=fliplr(FLIPFRAME);
            FLIPSF=fliplr(FLIPSF);
        end
        
        RDATA(:,:,frame)=FLIPFRAME;
        RSF=FLIPSF;
        waitbar(frame/size(NDATA,3))    
    end
    close(hdl);
else
    RDATA=NDATA;
    RSF=SF;
end    
%write final output results
DATA=RDATA;
SIGNAL=FR;
SF=RSF;
end    