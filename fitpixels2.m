%% DESCRIBE VARIABLES REQUIRED FROM MAIN PROGRAM
function [VFRAMEi,VFRAMEj,RMSEFRAME,XFRAMEINTERP_0,YFRAMEINTERP_0,SIGNALFRAMEINTERP_0]=fitpixels2(M,FrameXLim,FrameYLim,estvel,SIGNAL)
%% Define ROI
internal=1;%USE MATLAB2009 FIT FUNCTION
%estvel=mean(NEWDISTANCE)/dt;%[pixels/frame]
demopixel=0;%set to 1 to show fit
borderpixels=2;%extend frame area to provide neighbors for border pixels
pixelinterval=0.5;%interpolation
%cut frame if there is no room for borderpixel
if FrameXLim(1)<=borderpixels,FrameXLim(1)=3;end
if FrameXLim(2)>(size(M,2)-borderpixels),FrameXLim(2)=size(M,2)-borderpixels;end
if FrameYLim(1)<=borderpixels,FrameYLim(1)=3;end
if FrameYLim(2)>(size(M,1)-borderpixels),FrameYLim(2)=size(M,1)-borderpixels;end
    
XF=[FrameXLim(1)-borderpixels:FrameXLim(2)+borderpixels]';
YF=[FrameYLim(1)-borderpixels:FrameYLim(2)+borderpixels]';
[XFRAME,YFRAME]=meshgrid(XF,YF);
ZFRAME=M(YFRAME(1,1):YFRAME(end,1),XFRAME(1,1):XFRAME(1,end));%activation times

%Mark signals in frame+border
SIGNALFRAME=SIGNAL(YF(1):YF(end),XF(1):XF(end));%extended frame
%% interpolate activation times in frame
%create new grid
[XFRAMEINTERP,YFRAMEINTERP]=meshgrid([XF(1):pixelinterval:XF(end)],[YF(1):pixelinterval:YF(end)]);
%interpolate activation map
ZFRAMEINTERP=interp2(XFRAME,YFRAME,ZFRAME,XFRAMEINTERP,YFRAMEINTERP);
%interpolate SIGNALFRAME
SIGNALFRAMEINTERP=interp2(XFRAME,YFRAME,SIGNALFRAME,XFRAMEINTERP,YFRAMEINTERP);
SIGNALFRAMEINTERP(SIGNALFRAMEINTERP<1)=0;
%correlate SIGNALFRAMEINTERP with ZFRAMEINTERP
ZFRAMEINTERP(SIGNALFRAMEINTERP==0)=0;
%% fitting algorighm
%definitions outside look to save time

%initialize velocity components that include border
Vi=zeros(size(ZFRAMEINTERP));
Vj=zeros(size(ZFRAMEINTERP));
RMSE=zeros(size(ZFRAMEINTERP));

minpixels=6;%minimum number of pixels required for fit
minfitpixels=20;%minimum number of pixels desired for fit
%create fitoptions object
if internal==1
    OPTIONS = fitoptions('poly22');
end
%modify options structure
%lower boundaries
%L=[0,-Inf,-Inf,-Inf,-Inf,-Inf];
%U=[LIST(end),Inf,Inf,Inf,Inf,Inf];
% set(OPTIONS,'Lower',L,'Upper',U);
opts.Weights = zeros(1,0);

hdl = waitbar(0,'Fitting pixels');
for i=1:size(ZFRAMEINTERP,1)
    for j=1:size(ZFRAMEINTERP,2)
%% LOOP
        if SIGNALFRAMEINTERP(i,j)==1
        %time interval
        atime=ZFRAMEINTERP(i,j);
        %to improve wavefrontsearch, increase n by 1 if
        %length(NWFRONT)<minfitpixels
        repeatwavefrontsearch=1;n=1;
        while repeatwavefrontsearch==1 && n<=3
        
            SEARCHINTERVAL=[atime-n*minpixels/estvel/2,atime+n*minpixels/estvel/2];

            %find all pixels within the SEARCHINTERVAL
            WFRONT=[];
            for s=1:size(ZFRAMEINTERP,1)
                for t=1:size(ZFRAMEINTERP,2)
                    if SEARCHINTERVAL(1)<=ZFRAMEINTERP(s,t) && ZFRAMEINTERP(s,t)<=SEARCHINTERVAL(2)
                        if SIGNALFRAMEINTERP(s,t)==1
                            WFRONT=[WFRONT;[s,t,ZFRAMEINTERP(s,t)]];
                        end
                    end
                end
            end

        %define wavefront by reducing the number of WFRONT pixels to those that are
        %near the pixel of interest

            NWFRONT=[];
            for k=1:size(WFRONT,1)
                pixeldistance=norm([WFRONT(k,1)-i,WFRONT(k,2)-j],2);
                if pixeldistance<=minpixels/2+n-1
                    NWFRONT=[NWFRONT;WFRONT(k,:)];
                end
            end
        
        %if wavefront is too small (due to large increse in conduction velocity),
        %increase time interval
        
            if size(NWFRONT,1)<minfitpixels
                repeatwavefrontsearch=1;
                n=n+1;
            else
                repeatwavefrontsearch=0;
            end
        end
        
        %fit NWFRONT pixels
        %prepare fit
        %ALL POINTS
        X=NWFRONT(:,2);Y=NWFRONT(:,1);Z=NWFRONT(:,3);      
        %fit
        if size(NWFRONT,1)<minpixels
            %not enough points to perform fit
            RMSE(i,j)=0;
        else
            if internal==1
                %use fit function included in Matlab Ver>2009
                [fitresult, gof] = fit([X,Y],Z,'poly22',OPTIONS );
                RMSE(i,j)=gof.rmse;
            
                %calculate velocity components from fit
                %fitresult(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y +
                %p02*y^2
                p00=fitresult.p00;
                p10=fitresult.p10;
                p01=fitresult.p01;
                p20=fitresult.p20;
                p11=fitresult.p11;
                p02=fitresult.p02;

                %gradZ(x,y)=[Tx,Ty]
                Tx=p10+p11*i+2*p20*j;
                Ty=p01+p11*j+2*p02*i;
                Vx=Tx/(Tx^2+Ty^2);
                Vy=Ty/(Tx^2+Ty^2);
                Vi(i,j)=pixelinterval*Vy;%velocity in y-direction (row)
                Vj(i,j)=pixelinterval*Vx;%velocity in x-direction (column)
            else
                %use external fit function
                p=polyfitn([X(:),Y(:)],Z,2);
                %goodness of fit
                D=(Z(1)-polyvaln(p,[X(1),Y(1)]))^2;
                for n=2:length(Z)
                    D=D+(Z(n)-polyvaln(p,[X(n),Y(n)]))^2;
                end
                RMSE(i,j)=sqrt(D)/length(Z);
                P=p.Coefficients;
                %MODEL: T(x,y)=P(1)*X1^2 + P(2)*X1*X2 + P(3)*X1 + P(4)*X2^2 + P(5)*X2 + P(6)
                Tx=P(1)*2*j+P(2)*i+P(3);
                Ty=P(4)*2*i+P(2)*j+P(5);
                Vx=Tx/(Tx^2+Ty^2);
                Vy=Ty/(Tx^2+Ty^2);
                Vi(i,j)=pixelinterval*Vy;%velocity in y-direction (row)
                Vj(i,j)=pixelinterval*Vx;%velocity in x-direction (column)
            end
        end
    end
    end
    waitbar(i/size(ZFRAMEINTERP,1))
end
close(hdl)
%% prepare output of vector components
xmin=ceil(borderpixels/pixelinterval)+0.5;
ymin=xmin;
w=size(XFRAMEINTERP,2)-2*xmin+1-0.5;
h=size(YFRAMEINTERP,1)-2*ymin+1-0.5;
%crop velocity vector matrices           
VFRAMEi=imcrop(Vi,[xmin,ymin,w,h]);
VFRAMEj=imcrop(Vj,[xmin,ymin,w,h]);
RMSEFRAME=imcrop(RMSE,[xmin,ymin,w,h]);
XFRAMEINTERP_0=imcrop(XFRAMEINTERP,[xmin,ymin,w,h]);
YFRAMEINTERP_0=imcrop(YFRAMEINTERP,[xmin,ymin,w,h]);
SIGNALFRAMEINTERP_0=imcrop(SIGNALFRAMEINTERP,[xmin,ymin,w,h]);
framepixels=floor((bwarea(SIGNALFRAMEINTERP_0)));%total number of active pixels in frame
%% PLOT DATA FOR EXAMPLE PIXEL
if demopixel==1
%interpolate activation frame
PIXELij=[i,j];%PIXEL IN ZFRAMEINTERP (matrix coordinates)

atime=ZFRAMEINTERP(PIXELij(1),PIXELij(2));
f=10;%number of frames to consider for interpolation
df=0.1;%frame interval for interpolation
T=[floor(atime)-round(f/2):ceil(atime)+round(f/2)]';
Tinterp=[T(1):df:T(end)]';
i=1;while Tinterp(i)<=atime;i=i+1;end;aframe=i-1;%interpolated frame number

% %interpolate frames
% hdl = waitbar(0,'calculating activation frame');
% Ninterp=zeros(size(NORMDATA,1),size(NORMDATA,2),length(Tinterp));
% for m=1:size(NORMDATA,1)
%     for n=1:size(NORMDATA,2)
%         if SIGNAL(m,n)==1
%             P=squeeze(NORMDATA(m,n,T));
%             Ninterp(m,n,:)=interp1(T,P,Tinterp,'cubic');
%         end
%     end
%     waitbar(m/size(NORMDATA,1));
% end
% close(hdl);
% AFRAME=Ninterp(:,:,aframe);%corresponding frame with activation time closest to atime
% 
% %plot full depolarization map
% polfig=figure('name','POLARIZATION MAP');
% imagesc(AFRAME);hold on
% %mark frame + border
% rectangle('Position',[XFRAME(1,1)-0.5,YFRAME(1,1)-0.5,size(XFRAME,2),size(YFRAME,1)],'LineWidth',1,'EdgeColor','w','LineStyle','-');
% %mark original frame
% rectangle('Position',[XFRAME(1,1)-0.5+borderpixels,YFRAME(1,1)-0.5+borderpixels,size(XFRAME,2)-2*borderpixels,size(YFRAME,1)-2*borderpixels],'LineWidth',1,'EdgeColor','w','LineStyle','-');

% plot activation map
actfig=figure('name','ACTIVATION MAP');
imagesc(XFRAMEINTERP(1,:),YFRAMEINTERP(:,1),ZFRAMEINTERP);
rectangle('Position',[XFRAMEINTERP(1,PIXELij(2))-0.5*pixelinterval,YFRAMEINTERP(PIXELij(1),1)-0.5*pixelinterval,pixelinterval,pixelinterval],'LineWidth',1,'EdgeColor','w','LineStyle','-');

% plot wavefront
wavfig=figure('name','WAVEFRONT');
imagesc(ZFRAMEINTERP)
% plot wavefront pixels
for i=1:size(WFRONT,1)
    rectangle('Position',[WFRONT(i,2)-0.5,WFRONT(i,1)-0.5,1,1],'LineWidth',1,'EdgeColor','w','LineStyle','-')
end

% plot new wavefront
for i=1:size(NWFRONT,1)
    rectangle('Position',[NWFRONT(i,2)-0.5,NWFRONT(i,1)-0.5,1,1],'LineWidth',1,'EdgeColor','r','LineStyle','-')
end

% plot pixel of interest
rectangle('Position',[PIXELij(2)-0.5,PIXELij(1)-0.5,1,1],'LineWidth',1,'EdgeColor','k','LineStyle','-')

% plot fit
fitfig=figure('Name','FIT');
%plot DATA wavefront in BLACK
%X=NWFRONT(:,2);Y=NWFRONT(:,1);Z=NWFRONT(:,3);
plot3(X,Y,Z,'LineStyle','o','MarkerEdgeColor','k','MarkerFaceColor','k');

%plot fitted surface
[FITX,FITY]=meshgrid([min(X):0.2:max(X)],[min(Y):0.2:max(Y)]);  
FITZ=fitresult(FITX,FITY);
hold on
surf(FITX,FITY,FITZ,'EdgeColor','interp');

%plot fitted data points
if internal==1
    F=fitresult([X,Y]);
    fi=fitresult([PIXELij(2),PIXELij(1)]);
else
    F=polyvaln(p,[X,Y]);
    fi=polyvaln(p,[PIXELij(2),PIXELij(1)]);
end
plot3(X,Y,F,'LineStyle','o','MarkerEdgeColor','r','MarkerFaceColor','r');

%plot original point of interest
plot3(PIXELij(2),PIXELij(1),atime,'o','MarkerEdgeColor','m','MarkerFaceColor','m');
%plot fitted point
plot3(PIXELij(2),PIXELij(1),fi,'o','MarkerEdgeColor','m','MarkerFaceColor','m');

%plot fitted surface
% PZ=griddata(X,Y,Z,PX,PY);
end