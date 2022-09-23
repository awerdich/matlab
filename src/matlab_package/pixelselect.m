function [PIXELCOORDINATES,RGBFBR]=pixelselect(FBR,SCATTER,ASIGNALPIXELS,DATA)
%show image of the heart to select pixels
%convert signal image to RGB
%load bright field image 
if isempty(FBR)==1
[BFfile BFpath]=uigetfile('*.tif','Load image from high speed camera');
FBR=imread([BFpath,BFfile]);%load fluorescence image
end
FBR=imadjust(mat2gray(FBR));
[FBRIDX,map]=gray2ind(FBR,256);
RGBFBR=ind2rgb(FBRIDX,gray(256));

if isempty(ASIGNALPIXELS)==0 && isempty(SCATTER)==0
    SIGNALS=ASIGNALPIXELS;
    SIGNALS(SCATTER==0)=0;
    [SIGNALSIX,map]=gray2ind(SIGNALS,2);
    RGBSIGNALS=ind2rgb(SIGNALSIX,gray(2));
    ALPHASIG=SIGNALS*0.25;
    %change color of signals
    for i=1:size(SIGNALS,1)
        for j=1:size(SIGNALS,2)
            if SIGNALS(i,j)==1
                RGBSIGNALS(i,j,:)=[1,0,0];
            end
        end
    end
end

figure('name','Select two pixels at proximal and distal ends of heart');
%estimate pixels
%check if some pixels are already marked
image(RGBFBR);hold on
if isempty(ASIGNALPIXELS)==0 && isempty(SCATTER)==0
    image(RGBSIGNALS,'AlphaData',ALPHASIG)
    set(gca,'XLim',[1 80],'YLim',[1 80]);
end
[ax,ay]=ginput(2);
PIXELCOORDINATES=round([ax,ay]);