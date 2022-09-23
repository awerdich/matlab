%function [Vij,parameter,rmse,pdistance]=fitpixels(M,LOCATION,estvel,minfitpixels)
%% fitpixels function
displaycalculation=1;%show calculation for single pixel when run as script
scanrate=2000;
ERROR=[];PAR=[];DFRAME=[];DPIXEL=[];%store error and parameters for the three runs to decide which gives the better fit
i=LOCATION(1);
j=LOCATION(2);
nwfrontpixels=minfitpixels;%number of pixels to fit
rmse=0;%initialize rmse

%It is desirable that the ratio between the spatial (PIXELINTERVAL) and temporal window
%sizes (FRAMEINTERVAL) approximate propagation speed so that
%activity over the full spatial window will also span the temporal
%window
PIXELINTERVAL=[8,round(minfitpixels)+1];
dp=(PIXELINTERVAL(2)-PIXELINTERVAL(1))/10;
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
    SRMSE=[SRMSE;nrmse];
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