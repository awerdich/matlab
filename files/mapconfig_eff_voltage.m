%% configuration file for the optical mapping data analysis

%CAMERA CHIP DIMENSIONS CARDIO CCD-SMQ
pixelsize_x=24;%pixel size in um -> actual chip size divided by number of pixels 1920/80=24
pixelsize_y=24;

%SPECIAL CALIBRATIONS FOR C-MOUNT COUPLER COMBINATION 0.38X + 0.7X
%Calibration image H140730-10UMperDIV-O60S038&01-A.tif
%Nikon TE-2000
%Objective 60X Camera adapter 0.38 + 0.7 combined
%Camera SMQ
n=70;%number of pixels measured
d=100;%distance in um
obj=60;%objective magnification
p=d/n*60;%pixel-to-pixel distance without de-magnification
sideport1=pixelsize_x/p;

%SPECIAL CALIBRATIONS FOR C-MOUNT COUPLER COMBINATION 0.1X + 1.0X
%Camera Andor NEO
%with adapter 0.1 + 1.0
% n0=54;%number of pixels measured
% d0=100;%distance in um
% %without adapter 1.0 only
% n1=456;%number of pixels measured
% n2=467;%number of pixeels measured 2nd measurement
% d1=100;%distance in um
% sideport1=n0/mean([n1,n2]);%


%set default microscope configuration
if exist('magnification')==0 || exist('optovar')==0 || exist('sideport')==0 || exist('binning')==0
    magnification=40;
    optovar=1;
    sideport=0.28;
    binning=1;
end

%newmagnification=input(['Objective magnification [',num2str(magnification),']X >']);%objective magnification
%hard code to 20 for efficiency:
newmagnification=20

if isempty(newmagnification)==0 && newmagnification>0
    magnification=newmagnification;
end

%newoptovar=input(['Tube lens magnification (1/1.5) [',num2str(optovar),']X >']);%objective magnification
%hard code to 1:
newoptovar=1

if isempty(newoptovar)==0 && newoptovar>0
    if newoptovar==1.0 || newoptovar==1.5
        optovar=newoptovar;
    else
        error('Value not allowed (1.0 or 1.5 only).');
        clear magnification
    end
end

%fprintf(['0.38 + 0.7 sideport de-magnification ','= ',num2str(sideport1),'.\n']);
%newsideport=input(['Sideport magnification [',num2str(sideport),']X>']);%side port magnification
%hard code to 0.5 for efficiency:
newsideport=0.5

if isempty(newsideport)==0 && newsideport>0
    sideport=newsideport;
end

%newbinning=input(['Binning [',num2str(binning),'] >']);%camera hardware binning
%hard code to enter
newbinning=[]

if isempty(newbinning)==0 && newbinning>0
    binning=newbinning;
end

%pixel size
bin_x=binning;
bin_y=binning;
pixelcalfactor_x=pixelsize_x/(magnification*optovar*sideport)*bin_x;%final pixel size [um]
pixelcalfactor_j=pixelcalfactor_x;
pixelcalfactor_y=pixelsize_y/(magnification*optovar*sideport)*bin_y;%final pixel size [um]
pixelcalfactor_i=pixelcalfactor_y;

fprintf(['final interpixel distance (x,y): (',num2str(pixelcalfactor_x),',',num2str(pixelcalfactor_y),') um/pixel \n']);

%Screen size
P = get(0,'screensize');
screenwidth=P(1,3);
screenheight=P(1,4);
screensize=[screenwidth,screenheight];maxwindowsize=min(screensize(:));
windowsize=0.75.*[maxwindowsize,maxwindowsize];%[width,height]*screensize
offsetx=100;%distance between the left screen border and the first window
offsety=round((screenheight-windowsize(1))/2);%distance of windows from bottom of screen