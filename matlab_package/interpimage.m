%interpolate image
[imagefile imagepath]=uigetfile('*.tif','load image');
I=imread([imagepath,imagefile]);
J=imresize(I,10);
%save image
newfile=[imagefile(1:end-4),'-resize.tif'];
imwrite(J,[imagepath,newfile],'tif','Compression','none');