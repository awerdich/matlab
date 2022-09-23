%register background image to match fluorescence image
transparency=0.8;
%% load background image
if exist('backgroundfile')==0
    if exist('lastpath')==0
        lastpath=[];
    end
    [backgroundfile backgroundpath]=uigetfile('*.tif','load highres background image taken with observation camera',lastpath);    
    lastpath=backgroundpath;
    [backgroundfluorescencefile backgroundfluorescencepath]=uigetfile('*.tif','load highres fluorescence background image taken with observation camera or cancel',lastpath);
    if backgroundfluorescencefile~=0
        BF=imread([backgroundfluorescencepath backgroundfluorescencefile]);
        backgroundfluorescence=1;
    else
        BF=[];
        backgroundfluorescence=0;
    end
%% load bright field image for registration
    [BFfile BFpath]=uigetfile('*.tif','load brightfiled image taken with high speed camera',lastpath);
%% load fluorescence image
    [FLfile FLpath]=uigetfile('*.tif','load fluorescence image taken with high speed camera',lastpath);
    %FFL=squeeze(NORMDATA(:,:,1500));
end
%% load images
 FBR=imread([BFpath,BFfile]);%bright field image taken with fluorescence camera
 FFL=imread([FLpath,FLfile]);%fluorescence image taken with fluorescence camera
 B=imread([backgroundpath,backgroundfile]);
%% adjust ranges to [0 1] grayscale
GRAYB=mat2gray(B);
if isempty(BF)==0
    GRAYBF=mat2gray(BF);
else
    GRAYBF=GRAYB;
end
GRAYFBR=mat2gray(FBR);
GRAYFFL=mat2gray(FFL);
%% adjust contrast and filter background image
CB=imadjust(GRAYB);
CBF=imadjust(GRAYBF);

MCB=medfilt2(CB,[4 4]);
MCBF=medfilt2(CBF,[4 4]);
%% resize fluorescence and background images to 800 x 800 pixels

%resize fluorescence image
fresizefactor=10;%fluorescence image resize factor
RESFBR=imresize(GRAYFBR,fresizefactor,'bicubic');
RESFFL=imresize(GRAYFFL,fresizefactor,'bicubic');

%resize background image
zoom=1.3538565;%pixelsize ratio PIXELINK1024x768/REDSHIRT80x80
PIXELINKX=1024;%required imagesize in x direction for given zoom factor
if size(MCB,2)==PIXELINKX
    RB=imresize(MCB,zoom,'bicubic');
    RBF=imresize(MCBF,zoom,'bicubic');
else
%zoom factor does not match image size STOP PROGRAM
    fprintf(['wrong background image size for zoomfactor. PIXELNKX=',num2str(PIXELINKX),'!\n']);
    return
end

%display bright field image taken with fluorescence camera for comparison
figure('Name','Bright field image taken with fluorescence camera');
imagesc(RESFBR);colormap(gray);

figure('Name','select left upper edge of cropping area','Units','Pixels');
imagesc(RB);colormap(gray)

%set default cropping area
CENTERXY=[660,429];
ULXY=CENTERXY-400;
%show defaultcropping area
rectangle('Position',[ULXY(1),ULXY(2),800,800],'LineWidth',2,'EdgeColor','y');
rectangle('Position',[CENTERXY(1),CENTERXY(2),10,10],'LineWidth',2,'EdgeColor','y');

cropx=ULXY(1);cropy=ULXY(2);
if exist('cropx')==0
    [ax,ay,button]=ginput(1);
    if isempty(button)==1
        cropx=ULXY(1);cropy=ULXY(2);
    elseif button==1 || button==27
        cropx=round(ax);cropy=round(ay);
    end
end

%show cropping area
rectangle('Position',[cropx,cropy,800,800],'LineWidth',2,'EdgeColor','r');
rectangle('Position',[cropx+400,cropy+400,10,10],'LineWidth',2,'EdgeColor','r');

%crop image
RESB=imcrop(RB,[cropx,cropy,799,799]);
RESBF=imcrop(RBF,[cropx,cropy,799,799]);

%save frame coordinates
ulij=[cropy,cropx];%upper left corner in matrix coordinates
lrij=[cropy+size(RESB,1)-1,cropx+size(RESB,2)-1];%lower right corner in matrix coordinates


% [X,Y]=meshgrid(1:size(MCB,2),1:size(MCB,1));
% stepsize=size(MCB,1)/size(MCB,2);
% [XI,YI]=meshgrid(1:size(MCB,2),1:stepsize:size(MCB,1)+stepsize);
% RESB=interp2(X,Y,MCB,XI,YI,'cubic');
%% convert fluorescence image to rgb
[IFBR,mapFBR]=gray2ind(RESFBR,256);
[IFFL,mapFFL]=gray2ind(RESFFL,256);

RGBFBR=ind2rgb(IFBR,gray(256));
RGBFFL=ind2rgb(IFFL,gray(256));
%% image registration

if exist('LASTSHIFT')==1
    SHIFT=LASTSHIFT;
else 
    SHIFT=[0,0,0];
end

OFIG=figure('Name','IMAGE OVERLAY','MenuBar','none','Units','pixels');%prepare overlay figure
while isempty(SHIFT)==0

%shift background image RESB according to shift vector
BORDER=ones(size(RESB));%mark black border

if SHIFT(1)>=0 && SHIFT(2)>=0
    %horizontal shift of cropped backgound image
    if SHIFT(1)<cropx-1
        SHIFTBH=[RB(ulij(1):lrij(1),ulij(2)-SHIFT(1)+1:ulij(2)),RESB(:,1:end-SHIFT(1))];
        SHIFTBFH=[RBF(ulij(1):lrij(1),ulij(2)-SHIFT(1)+1:ulij(2)),RESBF(:,1:end-SHIFT(1))];
        BORDERH=BORDER;
    else
        SHIFTBH=[zeros(size(RESB,1),SHIFT(1)),RESB(:,1:end-SHIFT(1))];
        SHIFTBFH=[zeros(size(RESB,1),SHIFT(1)),RESBF(:,1:end-SHIFT(1))];
        BORDERH=[zeros(size(BORDER,1),SHIFT(1)),BORDER(:,1:end-SHIFT(1))];
    end
    %vertical shift of cropped background image
    if SHIFT(2)<size(RB,1)-lrij(1)
        SHIFTBV=[SHIFTBH(SHIFT(2)+1:end,:);RB(lrij(1):lrij(1)+SHIFT(2)-1,ulij(2)-SHIFT(1):lrij(2)-SHIFT(1))];
        SHIFTBFV=[SHIFTBFH(SHIFT(2)+1:end,:);RBF(lrij(1):lrij(1)+SHIFT(2)-1,ulij(2):lrij(2))];
        BORDERV=[BORDERH(SHIFT(2)+1:end,:);ones(SHIFT(2),size(BORDERH,2))];
    else
        SHIFTBV=[SHIFTBH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(SHIFTBH,2))];
        SHIFTBFV=[SHIFTBFH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(SHIFTBFH,2))];
        BORDERV=[BORDERH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(BORDERH,2))];
    end
    SHIFTB=SHIFTBV;
    SHIFTBF=SHIFTBFV;
    NBORDER=BORDERV;
elseif SHIFT(1)<0 && SHIFT(2)>=0
    %horizontal shift of cropped backgound image
    if abs(SHIFT(1))<size(RB,2)-lrij(2)
        SHIFTBH=[RESB(:,abs(SHIFT(1))+1:end),RB(ulij(1):lrij(1),lrij(2):lrij(2)+abs(SHIFT(1))-1)];
        SHIFTBFH=[RESBF(:,abs(SHIFT(1))+1:end),RBF(ulij(1):lrij(1),lrij(2):lrij(2)+abs(SHIFT(1))-1)];
        BORDERH=BORDER;
    else  
        SHIFTBH=[RESB(:,abs(SHIFT(1))+1:end),zeros(size(RESB,1),abs(SHIFT(1)))];
        SHIFTBFH=[RESBF(:,abs(SHIFT(1))+1:end),zeros(size(RESBF,1),abs(SHIFT(1)))];
        BORDERH=[BORDER(:,abs(SHIFT(1))+1:end),zeros(size(BORDER,1),abs(SHIFT(1)))];
    end
    %vertical shift of cropped background image
    if SHIFT(2)<size(RB,1)-(cropy+size(RESB,1)-1)
        SHIFTBV=[SHIFTBH(SHIFT(2)+1:end,:);RB(lrij(1):lrij(1)+SHIFT(2)-1,ulij(2)-SHIFT(1):lrij(2)-SHIFT(1))];
        SHIFTBFV=[SHIFTBFH(SHIFT(2)+1:end,:);RBF(lrij(1):lrij(1)+SHIFT(2)-1,ulij(2):lrij(2))];
        BORDERV=[BORDERH(SHIFT(2)+1:end,:);ones(SHIFT(2),size(BORDERH,2))];
    else
        SHIFTBV=[SHIFTBH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(SHIFTBH,2))];
        SHIFTBFV=[SHIFTBFH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(SHIFTBFH,2))];
        BORDERV=[BORDERH(SHIFT(2)+1:end,:);zeros(SHIFT(2),size(BORDERH,2))];
    end
    SHIFTB=SHIFTBV;
    SHIFTBF=SHIFTBFV;
    NBORDER=BORDERV;
elseif SHIFT(1)>=0 && SHIFT(2)<0
    %horizontal shift of cropped backgound image
    if SHIFT(1)<cropx-1
        SHIFTBH=[RB(ulij(1):lrij(1),ulij(2)-SHIFT(1)+1:ulij(2)),RESB(:,1:end-SHIFT(1))];
        SHIFTBFH=[RBF(ulij(1):lrij(1),ulij(2)-SHIFT(1)+1:ulij(2)),RESBF(:,1:end-SHIFT(1))];
        BORDERH=BORDER;
    else
        SHIFTBH=[zeros(size(RESB,1),SHIFT(1)),RESB(:,1:end-SHIFT(1))];
        SHIFTBFH=[zeros(size(RESB,1),SHIFT(1)),RESBF(:,1:end-SHIFT(1))];
        BORDERH=[zeros(size(BORDER,1),SHIFT(1)),BORDER(:,1:end-SHIFT(1))];
    end
    %vertical shift of cropped background image
    if abs(SHIFT(2))<ulij(1)-1
        SHIFTBV=[RB(ulij(1)-abs(SHIFT(2))+1:ulij(1),ulij(2):lrij(2));SHIFTBH(1:end-abs(SHIFT(2)),:)];
        SHIFTBFV=[RBF(ulij(1)-abs(SHIFT(2))+1:ulij(1),ulij(2):lrij(2));SHIFTBFH(1:end-abs(SHIFT(2)),:)];
        BORDERV=[ones(abs(SHIFT(2)),size(BORDERH,2));BORDERH(1:end-abs(SHIFT(2)),:)];
    else
        SHIFTBV=[zeros(abs(SHIFT(2)),size(SHIFTBH,2));SHIFTBH(1:end-abs(SHIFT(2)),:)];
        SHIFTBFV=[zeros(abs(SHIFT(2)),size(SHIFTBFH,2));SHIFTBFH(1:end-abs(SHIFT(2)),:)];
        BORDERV=[zeros(abs(SHIFT(2)),size(BORDERH,2));BORDERH(1:end-abs(SHIFT(2)),:)];
    end
    SHIFTB=SHIFTBV;
    SHIFTBF=SHIFTBFV;
    NBORDER=BORDERV;
elseif SHIFT(1)<0 && SHIFT(2)<0
   %horizontal shift of cropped backgound image
    if abs(SHIFT(1))<size(RB,2)-lrij(2)
        SHIFTBH=[RESB(:,abs(SHIFT(1))+1:end),RB(ulij(1):lrij(1),lrij(2):lrij(2)+abs(SHIFT(1))-1)];
        SHIFTBFH=[RESBF(:,abs(SHIFT(1))+1:end),RBF(ulij(1):lrij(1),lrij(2):lrij(2)+abs(SHIFT(1))-1)];
        BORDERH=BORDER;  
    else  
        SHIFTBH=[RESB(:,abs(SHIFT(1))+1:end),zeros(size(RESB,1),abs(SHIFT(1)))];
        SHIFTBFH=[RESBF(:,abs(SHIFT(1))+1:end),zeros(size(RESBF,1),abs(SHIFT(1)))];
        BORDERH=[BORDER(:,abs(SHIFT(1))+1:end),zeros(size(BORDER,1),abs(SHIFT(1)))];
    end
    %vertical shift of cropped background image
    if abs(SHIFT(2))<ulij(1)-1
        SHIFTBV=[RB(ulij(1)-abs(SHIFT(2))+1:ulij(1),ulij(2):lrij(2));SHIFTBH(1:end-abs(SHIFT(2)),:)];
        SHIFTBFV=[RBF(ulij(1)-abs(SHIFT(2))+1:ulij(1),ulij(2):lrij(2));SHIFTBFH(1:end-abs(SHIFT(2)),:)];
        BORDERV=[ones(abs(SHIFT(2)),size(BORDERH,2));BORDERH(1:end-abs(SHIFT(2)),:)];
    else
        SHIFTBV=[zeros(abs(SHIFT(2)),size(SHIFTBH,2));SHIFTBH(1:end-abs(SHIFT(2)),:)];
        SHIFTBFV=[zeros(abs(SHIFT(2)),size(SHIFTBFH,2));SHIFTBFH(1:end-abs(SHIFT(2)),:)];
        BORDERV=[zeros(abs(SHIFT(2)),size(BORDERH,2));BORDERH(1:end-abs(SHIFT(2)),:)];
    end
    SHIFTB=SHIFTBV;
    SHIFTBF=SHIFTBFV;
    NBORDER=BORDERV;
end

%rotate image according to SHIFT(1,3)
if length(SHIFT)<3
    SHIFT=[SHIFT,0];
end
SHIFTB=imrotate(SHIFTB,SHIFT(3),'bicubic','crop');
SHIFTBF=imrotate(SHIFTBF,SHIFT(3),'bicubic','crop');
NBORDER=imrotate(NBORDER,SHIFT(3),'crop');   

%convert background image to rgb
[IB,mapB]=gray2ind(SHIFTB,256);%convert 16 bit intensity image into 16 bit indexed image 
RGBSHIFTB=ind2rgb(IB,hot(256));

[IBF,mapBF]=gray2ind(SHIFTBF,256);
RGBSHIFTBF=ind2rgb(IB,jet(256));

%display both images
%background
figure(OFIG);
image(RGBSHIFTB);
hold on;
%display bright field image taken with high speed camera
ALPHABR=alphamatrix(RESFBR,0,transparency);
image(RGBFBR,'AlphaData',ALPHABR);
hold off

%ask for shift parameters
LASTSHIFT=SHIFT;
SHIFT=input(['pixels to shift [x,y,angle]: [',num2str(SHIFT(1)),',',num2str(SHIFT(2)),',',num2str(SHIFT(3)),'] :']);
SHIFT=round(SHIFT);
end
SHIFT=LASTSHIFT;
%% plot shifted background image with fluorescence overlay
figure;
image(RGBSHIFTB);hold on
ALPHAFL=ones(size(RESFFL))*0.6;
%ALPHAFL=alphamatrix(RESFFL,0.2,0.5);
image(RGBFFL,'AlphaData',ALPHAFL);
%% color black areas resulting from image registration
figure('name','click in background');imagesc(SHIFTB);colormap(gray)
[bx,by,button]=ginput(1);

%remove decimals from bx and by
ctri=sprintf('%2.1f',by);ctrj=sprintf('%2.1f',bx);
ctri=str2num(ctri(1:end-2));ctrj=str2num(ctrj(1:end-2)); 
meanctr=mean(SHIFTB(ctri-5:ctri+5,ctrj));

SHIFTB(NBORDER==0)=meanctr; 
SHIFTBF(NBORDER==0)=meanctr;
imagesc(SHIFTB);colormap(gray);

%rescale background images to [0 1]
SHIFTB=mat2gray(SHIFTB);
SHIFTBF=mat2gray(SHIFTBF);

%re-set NBORDER pixels in background fluorescence image
SHIFTBF(NBORDER==0)=0;
%convert new background image to rgb
[IB,mapB]=gray2ind(SHIFTB,256);%convert 16 bit intensity image into 16 bit indexed image 
RGBSHIFTB=ind2rgb(IB,gray(256));
[IBF,mapBF]=gray2ind(SHIFTBF,256);
RGBSHIFTBF=ind2rgb(IBF,gray(256));

%assign color to background fluorescence
RGBSHIFTBF(:,:,[1,3])=0;%take only green color component

%% save registered image
BRfile=[backgroundfile(1:end-4),'-RegBackground.mat'];
if exist('datapath')==0
    [file datapath]=uiputfile('*.*','SELECT FOLDER');
end
if isempty(BF)==1
    SHIFTBF=[];
    RGBSHIFTBF=[];
end 
save([datapath,BRfile],'B','BF','FBR','FFL','NBORDER','meanctr','SHIFTB','SHIFTBF','RGBSHIFTB','RGBSHIFTBF','SHIFT','fresizefactor')