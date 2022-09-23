%program to overlay camera pixel array with brightfield image
[cfile cpath]=uigetfile('*.tif','OPEN CAMER IMAGE FILE');
I=imread([cpath,cfile]);
arrayscalefactor=200/61;%um/pixels
camerascalefactor=220/743;%um/pixels
arraypixel=arrayscalefactor/camerascalefactor;%array pixel size in camera pixels
arrayelements=80;%number of array elements
cameraelements=[800,600];%camera resolution

fig = figure('Name','Bright field image of the heart','MenuBar','none','Units','pixels','Position',[300 400 800 600],'Color','k','Visible','on');50
heart=axes('Units','normalized','Position',[0 0 1 1],'Visible','on','Color','w');

axes(heart),set(heart,'TickDir','out','Units','pixels'),image(I);

%add scalebar to image
scalebarlength=50;%length of scalebar in picture (um)
scalebarpixels=scalebarlength/camerascalefactor;%length of scalebar (pixels)
unitstring=[num2str(scalebarlength),' \mum'];

scalebarxdata=(cameraelements(1)-scalebarpixels):cameraelements(1);
scalebarydata=cameraelements(2)-10*ones(1,length(scalebarxdata));

%draw scalebars at the bottom right corner of image
line('xdata',scalebarxdata,'ydata',scalebarydata,'color','w','LineWidth',2);

%add unit length to horizontal scalebar
text(cameraelements(1)-110,cameraelements(2)-25,unitstring,'color','w','FontWeight','bold')


%select pixel position
[xi,yi]=ginput(1);
P=[round(xi)-arraypixel/2,round(yi)-arraypixel/2];

axes(heart),rectangle('Position',[P(1),P(2),arraypixel,arraypixel],'EdgeColor','r','LineWidth',2);


%save image
%get current image and save as movie frame
F=getframe(gcf);
%convert movieframe back into image
[IMAGE,imagemap]=frame2im(F);
%write image to file
imagefile=[cfile(1:end-4),'.jpg'];
imagepath=cpath;
imwrite(IMAGE,[imagepath,imagefile],'jpg','Quality',100);
