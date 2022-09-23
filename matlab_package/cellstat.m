%compute statistics of isolated objects in an image
if exist('setumlength')==0 || isempty(setumlength)==1
    setumlength=210.42;
end
if exist('setpxlength')==0 || isempty(setpxlength)==1
    setpxlength=512;
end
if exist('setrefangle')==0 || isempty(setrefangle)==1
    setrefangle=180;
end

pxlength=input(['image resolution (pixels) [',num2str(setpxlength),'] >']);
umlength=input(['calibrated length (microns) [',num2str(setumlength),'] >']);
refangle=input(['reference angle (degrees) [',num2str(setrefangle),'] 180=enter from image>']);

if isempty(pxlength)==1
    pxlength=setpxlength;
else
    setpxlength=pxlength;
end
if isempty(umlength)==1
    umlength=setumlength;
else
    setumlength=umlength;
end
if isempty(refangle)==1
    refangle=setrefangle;
else
    setrefangle=refangle;
end

pixelcalfactor=umlength/pxlength;
fprintf(['pixel calibration: ',num2str(pixelcalfactor),' um/pixel\n']);
%% read and display image
if exist('imagefile')==0 || isempty(imagefile)==1
    [imagefile imagepath]=uigetfile('*.tif','LOAD IMAGE');
end
I1=imread([imagepath,imagefile]);%load image
I2=imadjust(I1);%improve contrast
figure;imshow(I2);hold on
%enter reference angle manually if refangle=180
if refangle==180
    fprintf('Draw reference line in image. [RETURN] to proceed. \n')
    h = imline(gca,[],[]);
    api = iptgetapi(h);
    fcn = makeConstrainToRectFcn('imline',get(gca,'XLim'),get(gca,'YLim'));
    api.setPositionConstraintFcn(fcn);
    id = api.addNewPositionCallback(@(pos) title(mat2str(pos,3)));
    %precent dragging line outside extent of image
    pause
    api.removeNewPositionCallback(id);
    position=api.getPosition();
    ax=position(:,1);ay=position(:,2);
    line(ax,ay,'color','r','LineWidth',1.5);
    line(ax,[max(ay);max(ay)],'color','k','LineWidth',1.5);
    R=[ax(2);ay(2)]-[ax(1);ay(1)];
    if ay(1)<ay(2)
        R=-R;
    end
    refangle=acos(R(1)/norm(R))*180/pi;
    text(min(ax),max(ay)-10,['angle:',num2str(round(refangle)),'deg']);
    fprintf(['Reference angle: ',num2str(refangle),' deg. [RETURN] to proceed\n']);
    pause
end
%% crop and convert grayscale image to binary binary image
%crop image
fprintf('Crop image. Draw rectangle in image and select crop from menu.\n');
I3=imcrop(I2);
%adjust range to [0,1] 
I4=mat2gray(I3);
%convert grayscale [0 1] image to binary image using global thresholding
level=graythresh(I4);
I5=im2bw(I4,level);
%invert bw image so that the objects are 1
I6=~I5;
%fill objects in binary image
I7=imfill(I6,'holes');
%% define objects in image
minpixel=10;%minimum number of pixels that make up an object
[LABEL,num]=bwlabel(I7,8);
I8=I7;%modify copy of original image
%cound number of pixels in each object and delete all objects in I8 that
%are smaller than minpixel

DELOBJ=[];%[o-number,pixels] number of objects removed from image
for k=1:num
    %get pixel coordinates for object i
    [I,J] = find(LABEL==k);
    %delete object if it is too small
    if length(I)<minpixel
        for m=1:length(I)
            I8(I(m),J(m))=0;
        end
        DELOBJ=[DELOBJ;[k,length(I)]];%number of deleted object
    end
end
%re-count objects
[NEWLABEL,newnum]=bwlabel(I8,8);
%% measure object properties
REGION=NEWLABEL;%LABEL MATRIX
CELLSTAT=regionprops(REGION,'all');
%% calculate longest difference vector for the eight extrema points
for k=1:newnum
    E=CELLSTAT(k).Extrema;
    VECTORLIST=[];
    for i=1:length(E)-1
        for j=(i+1):length(E)
            D=[E(j,1);E(j,2)]-[E(i,1);E(i,2)];
            %invert vector if it points downwards;
            if D(2)>0
                D=-D;
            end 
            angle=acos(D(1)/norm(D))*180/pi;
            if norm(D)>0
                VECTORLIST=[VECTORLIST;[D(1),D(2),angle,norm(D)]];
            end
        end
    end
    %sort vectorlist to find longest axis
    [S,SX]=sort(VECTORLIST(:,4),'ascend');
    %add longest vector and angle to CELLSTAT
    CELLSTAT(k).Axis=[VECTORLIST(SX(end),1),VECTORLIST(SX(end),2)];
    CELLSTAT(k).Angle=VECTORLIST(SX(end),3);
end
%% calculate circularity
for k=1:newnum
    CELLSTAT(k).Circularity=4*pi*CELLSTAT(k).Area/CELLSTAT(k).Perimeter^2;
end 
%% prepare output image
%convert binary image to RGB image
[I8IX,I8MAP]=gray2ind(I8,256);
I8RGB=ind2rgb(I8IX,jet(256));
%color background and objects
[Bi,Bj]=find(REGION==0);%coordinates of background
[Oi,Oj]=find(REGION>0);%coordinates of all objects
%color background
for i=1:length(Bi)
    I8RGB(Bi(i),Bj(i),:)=[1,1,1];%WHITE
end
%color objects
for i=1:length(Oi)
    I8RGB(Oi(i),Oj(i),:)=[1,0,0];%RED
end
%% draw measurements
figure;image(I8RGB);hold on
for k=1:length(CELLSTAT)
    %center of mass
    COM=[CELLSTAT(k).Centroid(1),CELLSTAT(k).Centroid(2)];
    %draw center of mass
    rectangle('Position',[round(COM(1))-0.5,round(COM(2))-0.5,1,1],'FaceColor','b')
    %annotate object
    text(round(COM(1))+1.5,round(COM(2))+1.5,num2str(k))
    %vector 
    vx=CELLSTAT(k).Axis(1);vy=CELLSTAT(k).Axis(2);
    quiver(COM(1),COM(2),vx,vy,0.5,'Color','b','LineWidth',1.5,'MaxHeadSize',1.5)
end
%% save data
datafile=[imagefile(1:end-4),'.mat'];
save([imagepath datafile],'CELLSTAT','I1','I8','I8RGB','umlength','pxlength','pixelcalfactor','imagefile','refangle');
%prepare ascii data
CELLDATA=[];
%CELLDATA=[CELL,AREA,PERIMETER,CIRCULARITY,ANGLE]
for k=1:length(CELLSTAT)
    kangle=CELLSTAT(k).Angle-refangle;
    if kangle<0
        kangle=180+kangle;
    end
    S=[k,CELLSTAT(k).Area*pixelcalfactor^2,CELLSTAT(k).Perimeter*pixelcalfactor,CELLSTAT(k).Circularity,kangle];
    CELLDATA=[CELLDATA;S];
end
datatextfile=[imagefile(1:end-4),'.txt'];
save([imagepath datatextfile],'CELLDATA','-ascii');