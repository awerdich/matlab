%% load background image
if exist('backgroundfile')==0
    [backgroundfile backgroundpath]=uigetfile('*.tif','load highres background image taken with observation camera');    
    B=imread([backgroundpath,backgroundfile]);
    [backgroundfluorescencefile backgroundfluorescencepath]=uigetfile('*.tif','load highres fluorescence background image taken with observation camera or cancel');
    BF=imread([backgroundfluorescencepath backgroundfluorescencefile]);   
end
%% re-scale pixel values
GRAYB=mat2gray(B);
GRAYBF=mat2gray(BF);
%% adjust contrast and filter image
CB=imadjust(GRAYB);
CBF=imadjust(GRAYBF);

MCB=medfilt2(CB,[4 4]);
MCBF=medfilt2(CBF,[4 4]);
%% convert fluorescence image to rgb
[IB,mapIB]=gray2ind(MCB,256);
[IBF,mapIBF]=gray2ind(MCBF,256);

RGBB=ind2rgb(IB,gray(256));
RGBF=ind2rgb(IBF,gray(256));
% change color of fluorescence image
for i=1:size(RGBF,1)
    for j=1:size(RGBF,2)
        RGBF(i,j,1)=0;
        RGBF(i,j,3)=0;
    end
end
%% create overlay image
% figure and axes
exportfigure = figure('Name','HEART VIEW','MenuBar','none','Units','pixels','Position',[100 10 1024 768],'Color','k','Visible','on');
exportaxes=axes('Position',[0 0 1 1],'Visible','on','Drawmode','fast');
% Alphamatrix
ALPHA=alphamatrix(MCBF,0,0.6);

% display images
axes(exportaxes);
image(RGBB);
hold on
image(RGBF,'AlphaData',ALPHA);
set(exportaxes,'TickDir','Out');
%% save frame as tif file
%get current image and save as movie frame
figure(exportfigure);
F=getframe(gcf);
%convert movieframe back into image
[IMAGEC,imagemap]=frame2im(F);
%crop image to remove black border
%IMAGEC=imcrop(IMAGE,[2,2,size(IMAGE,1)-4,size(IMAGE,2)-4]);
%write image to file
[imagefile,imagepath]=uiputfile('*.tif','pick path');
region=input('description image>','s');
imagefile=[stackfile(1:end-4),'-',region,'.tif'];
imagepathfile=[imagepath,imagefile];
imwrite(IMAGEC,imagepathfile,'tif','Compression','none');
