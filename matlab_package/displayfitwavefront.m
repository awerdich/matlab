%% parameters
M=C;
%combine ASIGNALPIXELS and SCATTER
NSIGNALPIXELS=ASIGNALPIXELS;
if exist('SCATTER')==1
    for k=1:size(SCATTER,1)
        NSIGNALPIXELS(SCATTER(k,1),SCATTER(k,2))=0;
    end
end
LOCATION=[i,j];
minfitpixels=10;%minimum number of pixels used in fit
estvelmms=15.0;%estimated propagation speed [mm/s]
%% fitpixels function
scanrate=2000;
pixelcalfactor=16/7.201613;%REDSHIRT 80x80 [um/pixel]
calvfactor=scanrate*1.0e-3*pixelcalfactor;%[frame*mm/(pixel*s)][pixel/frame*frame/s*um/pixel*1000mm/um=mm/s]
estvel=estvelmms/calvfactor;
ERROR=[];PAR=[];DFRAME=[];DPIXEL=[];%store error and parameters for the three runs to decide which gives the better fit
i=LOCATION(1);
j=LOCATION(2);
nwfrontpixels=minfitpixels;%number of pixels to fit
rmse=0;%initialize rmse

%It is desirable that the ratio between the spatial (PIXELINTERVAL) and temporal window
%sizes (FRAMEINTERVAL) approximate propagation speed so that
%activity over the full spatial window will also span the temporal
%window
PIXELINTERVAL=[2,10];
dp=0.5;
%initialize results
SPARAMETER=[];
SRMSE=[];
SVij=[];
SPIXELDISTANCE=[];
SWFRONT=[];
SNWFRONT=[];
SWFRONTSIZE=[];
SNWFRONTSIZE=[];
for tc=1:3
    maximumpixeldistance=PIXELINTERVAL(1);
while maximumpixeldistance<PIXELINTERVAL(2)
    dt=maximumpixeldistance/estvel;
    %initialize variables
    %save wavefront from previous iteration
    
    WFRONT=[];%pixels with activation times within the searchinterval
    NWFRONT=[];%pixels in WFRONT that are minimumpixeldistance pixels away from the fitted pixel
    %find neighbors to (i,j) on the wavefront
    atime=M(i,j);%activation time of pixel (i,j)  
    %USE DIFFERENT SEARCHINTERVALS (all dt wide) FOR 3 RUNS tc
    if tc==1
        SEARCHINTERVAL=[atime-dt,atime];%corresponding time interval in frames
    elseif tc==2
        SEARCHINTERVAL=[atime-dt/2,atime+dt/2];
    elseif tc==3
        SEARCHINTERVAL=[atime,atime+dt];
    end
    %build wavefront
    WFRONT=[];
    for t=1:size(M,2)
        for s=1:size(M,1)
            if SEARCHINTERVAL(1)<=M(s,t) && M(s,t)<=SEARCHINTERVAL(2)
                WFRONT=[WFRONT;[s,t,M(s,t)]];
            end
        end
    end
    %throw away pixels on the wavefront that are too far away
    
    NWFRONT=[];
    for k=1:size(WFRONT,1)
        DIFFVECTOR=[WFRONT(k,1)-i;WFRONT(k,2)-j];
        pixeldistance=norm(DIFFVECTOR,2);%2-norm of DIFFVECTOR
        if pixeldistance<maximumpixeldistance
            NWFRONT=[NWFRONT;[WFRONT(k,1),WFRONT(k,2),WFRONT(k,3),k]];
        end
    end
if size(NWFRONT,1)<nwfrontpixels
    %if not enough pixels found, try to increase maximum pixeldistance
    maximumpixeldistance=maximumpixeldistance+dp;
    success=0;
else
   %required number of pixels on wavefront found
   savedistance=maximumpixeldistance;
   maximumpixeldistance=PIXELINTERVAL(2); 
   success=1;
end
%end

% FIT
if success==1
X=NWFRONT(:,2);Y=NWFRONT(:,1);Z=NWFRONT(:,3);

    %save parameters from previous successful fit to compare with new fit
    
    
    %hyperbolic fit
    p=polyfitn([X(:),Y(:)],Z,2);
    
    %estimate goodness of fit
    D=(Z(1)-polyvaln(p,[X(1),Y(1)]))^2;
    for n=2:length(Z)
        D=D+(Z(n)-polyvaln(p,[X(n),Y(n)]))^2;
    end
    nrmse=sqrt(D)/length(Z);

    %if this fit was worse than the last one, keep last one
          
elseif success==0
    nrmse=0;
end
end

% calculate velocity components using the fit
if nrmse>0
    x=LOCATION(2);y=LOCATION(1);P=p.Coefficients;
    %MODEL: T(x,y)=P(1)*X1^2 + P(2)*X1*X2 + P(3)*X1 + P(4)*X2^2 + P(5)*X2 + P(6)
    Tx=P(1)*2*x+P(2)*y+P(3);
    Ty=P(4)*2*y+P(2)*x+P(5);
    Vx=Tx/(Tx^2+Ty^2);
    Vy=Ty/(Tx^2+Ty^2);
    Vij(1)=Vy;%velocity in y-direction (row)
    Vij(2)=Vx;%velocity in x-direction (column)
else
    Vij=[0,0];
end
%save fit results for each run, if fit was successul
if nrmse>0
    SPARAMETER=[SPARAMETER;p];
    SRMSE=[SRMSE;nrmse]
    SVij=[SVij;Vij];
    SPIXELDISTANCE=[SPIXELDISTANCE;maximumpixeldistance];
    SWFRONT=[SWFRONT;WFRONT];
    SNWFRONT=[SNWFRONT;NWFRONT];
    SWFRONTSIZE=[SWFRONTSIZE;size(WFRONT,1)];
    SNWFRONTSIZE=[SNWFRONTSIZE;size(NWFRONT,1)];
end
end

%pick best result if at least one fit was successful
if isempty(SRMSE)==0
    [val,idx]=min(SRMSE);
    parameter=SPARAMETER(idx);
    rmse=SRMSE(idx);
    Vij=SVij(idx,:);
    pdistance=SPIXELDISTANCE(idx);
    if idx>1
        BESTWFRONT=SWFRONT(sum(SWFRONTSIZE(1:idx-1))+1:sum(SWFRONTSIZE(1:idx-1))+SWFRONTSIZE(idx),:);
        BESTNWFRONT=SNWFRONT(sum(SNWFRONTSIZE(1:idx-1))+1:sum(SNWFRONTSIZE(1:idx-1))+SNWFRONTSIZE(idx),:);
    else
        BESTWFRONT=SWFRONT(1:SWFRONTSIZE(idx),:);
        BESTNWFRONT=SNWFRONT(1:SNWFRONTSIZE(idx),:);
    end
else %if none of the fits was successful return rmse=0 
    parameter=[];
    rmse=0;
    Vij=[];
    pdistance=[];
end

%% display frame and pixels on wavefront corresponding to (i,j) in POL
WFRONT=BESTWFRONT;
NWFRONT=BESTNWFRONT;
%define axes for picture display
heartfigure = figure('Name','HEART VIEW','Units','normalized','Position',[0.25 0.2 0.5 0.55],'Color','w','Visible','on');
%set axes for figheart
figure(heartfigure);
heart=axes('Position',[0.1 0.1 0.8 0.8],'Visible','on','Drawmode','fast');

%round activation time to closest frame and display frame
frame=round(M(i,j));
I=NORMDATA(:,:,frame);
I(NSIGNALPIXELS==0)=0;
%convert grayscale to rgb image
[XI,map]=gray2ind(I,256);%convert 16 bit intensity image into 16 bit indexed image 
J=ind2rgb(XI,jet(256));%convert indexed image into truecolor image using colormap jet
%change background to white
for k=1:size(NSIGNALPIXELS,1)
    for l=1:size(NSIGNALPIXELS,2)
        if NSIGNALPIXELS(k,l)==0
            J(k,l,:)=[1 1 1];
        end
    end
end

%plot frame
axes(heart),image(J)
set(heart,'TickDir','out','XLim',[1 80],'YLim',[1 80])

%plot wavefront pixels
for n=1:size(WFRONT,1)
    rectangle('Position',[WFRONT(n,2)-0.5,WFRONT(n,1)-0.5,1,1],'LineWidth',1.5,'EdgeColor','w');
end

%plot wavefront pixels that are in the neighborhood
for n=1:size(NWFRONT,1)
    rectangle('Position',[NWFRONT(n,2)-0.5,NWFRONT(n,1)-0.5,1,1],'LineWidth',1.5,'EdgeColor','r');
end

%plot active pixel
rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',1.5,'EdgeColor','k');

%% display fit of wavefront to surface
%display fit result
%generate data on regular grid
%FITZ=feval(@wavefit,P,DATA);

%generate uniformly spaced data with the desired number of elements
XLIN=linspace(min(X),max(X),100);
YLIN=linspace(min(Y),max(Y),100);
%using these points, generate a uniformly spaced grid
[IX,IY]=meshgrid(XLIN,YLIN);

%evaluate fitted surface polynom at the uniformly spaced points
IZ=zeros(size(IX));

for m=1:size(IX,2)
    for n=1:size(IY,1)
        a=IX(1,m);b=IY(n,1);
        IZ(n,m)=polyvaln(parameter,[a b]);
    end
end

%plot fitted surface and data points
figure;
plot3(X,Y,Z,'b.','MarkerSize',20);hold on;
%plot3(X,Y,polyvaln(parameter,[X,Y]),'b+','MarkerSize',20);
plot3(j,i,M(i,j),'r.','MarkerSize',20);
plot3(j,i,polyvaln(parameter,[j,i]),'r.','MarkerSize',20);
axis tight;
surf(IX,IY,IZ)