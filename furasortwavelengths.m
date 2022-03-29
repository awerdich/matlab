function [W1,W2]=furasortwavelengths(DATAST)
%DATAPATH=stackpath;DATAFILE=stackfile;
%open data file frame-by-frame
%get total number of frames
width=DATAST(1).columns;
height=DATAST(1).rows;
numframes=size(DATAST(1).RAWDATA,3);
RAWDATA=DATAST(1).RAWDATA;
fratio=DATAST(1).fratio;
%% Sort Wavelengths
maxframe=floor(numframes/fratio)*fratio;%highest frame number for complete ratios
numratios=maxframe/fratio - 1;%number of complete ratios
W1=zeros(size(RAWDATA,1),size(RAWDATA,2),numratios);%first wavelength intensities
W1FRAME=zeros(numratios,1);%absolute frame numbers of images in W1
W2=zeros(size(W1));%second wavelength intensities
W2FRAME=zeros(size(W1FRAME));
f1 = [];
f2 = [];
offset = 0;
for ratio=1:numratios
    frame1=fratio/2*(2*ratio-1) + offset;%
    f1 = [f1; frame1];
    W1(:,:,ratio)=RAWDATA(:,:,frame1);
    W1FRAME(ratio)=frame1;
    frame2=fratio*ratio + offset;%EVEN FRAMES
    f2 = [f2; frame2];
    W2(:,:,ratio)=RAWDATA(:,:,frame2);
    W2FRAME(ratio)=frame2;
end
