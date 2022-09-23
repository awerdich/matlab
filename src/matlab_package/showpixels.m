%Plot activated pixels into potential figure

DATA=NORMDATA;
ACTIVATIONMATRIX=POL;
%ACTIVATIONMATRIX=UPSTROKEIDX;

%choose background color
%bcolor=input('background color (w) white (k) black (r) red (g) green (b) blue >','s');
bcolor='w';


if bcolor=='w'
   brgb=[1 1 1];
   tx='k';
elseif bcolor=='k'
    brgb=[0 0 0];
    tx='w';
elseif bcolor=='r'
    brgb=[1 0 0];
    tx='w';
elseif bcolor=='g'
    brgb=[0 1 0];
    tx='k';
elseif bcolor=='b'
    brgb=[0 0 1];
    tx='w';
end
    


%define axes for picture display

heartfigure = figure('Name','HEART VIEW','MenuBar','none','Units','normalized','Position',[0.25 0.41 0.5 0.55],'Color','w','Visible','on');

%set axes for figheart
figure(heartfigure);
heart=axes('Position',[0.1 0.1 0.8 0.8],'Visible','on','Drawmode','fast');

%display frames

frame=1;
button=0;
colormap jet(65535);

%choose frame
while isempty(frame)==0    
    
    %plot frame
    I=DATA(:,:,frame);
    
    %find background in frame
    B=[];
    for i=1:size(I,1)
        for j=1:size(I,2)
            if mean(DATA(i,j,frame:frame+10))==0
                B=[B;[i,j]];
            end
        end
    end
    
    [XI,map]=gray2ind(I,256);%convert 16 bit intensity image into 16 bit indexed image 
    J=ind2rgb(XI,jet(256));%convert indexed image into truecolor image using colormap jet
    %stain background
    for i=1:size(B,1),J(B(i,1),B(i,2),:)=brgb;,end

%% Plot heart
   
    %Plot heart
    axes(heart),image(J),set(heart,'TickDir','out')

%% search pixels in activationmatrix by given interval
peaknumber=1;
cframe=frame+0;
SEARCHINTERVAL=[cframe-1 cframe+1];
C=[];
for i=1:size(ACTIVATIONMATRIX,1)
    for j=1:size(ACTIVATIONMATRIX,2)
        if SEARCHINTERVAL(1)<=ACTIVATIONMATRIX(i,j,peaknumber) && ACTIVATIONMATRIX(i,j,peaknumber)<=SEARCHINTERVAL(2)
            C=[C;[i,j]];
        end
    end
end
  
%add frames for pixels in matrix C
axes(heart)
for i=1:size(C,1)
    rectangle('Position',[C(i,2)-0.5,C(i,1)-0.5,1,1],'LineWidth',1.5,'EdgeColor','w');
end
    
%% imput new frame

    %input frame number
    oldframe=frame;
    frame=input(['Current frame: ',num2str(frame),'/',num2str(size(DATA,3)),' frame>']);
end

frame=oldframe;%remember last frame