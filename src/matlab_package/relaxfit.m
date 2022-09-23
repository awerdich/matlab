%% relaxation fit to sarcomere data
%data
CTRANS=PL_F;%calcium transients
STRANS=SARC_PL;%sarcomere contraction transients
REST=SARC_PL_REST;%resting sarcomere length
stimulus_scan=600;%common stimulus upstroke scan
BASELINE=[1:500]';
scanrate=2.0e4;%scan rate (Hz)
%% fit base line to compensate for photo bleaching
%fit baseline
FITCA=zeros(size(CTRANS));
XDATA=[BASELINE;[6000:8000]'];
for i=1:size(CTRANS,2)
    CA=CTRANS(:,i);
    P=polyfit(XDATA,CA(XDATA),1);
    FITCA(:,i)=CA-((polyval(P,1:length(CA))')-mean(BASELINE));
end
%% correct offset and normalize 
NCA=zeros(size(FITCA));
for i=1:size(FITCA,2)
    CA=FITCA(:,i)-mean(FITCA(BASELINE,i),1);
    NCA(:,i)=CA/max(CA);
end
%% find 95% ca uptake
threshold=0.04;
FITLEVELS=[0.4,0.1];
INDEX95=zeros(size(NCA,2),1);
for i=1:size(NCA,2)
    T=NCA(:,i);
    %determine maximum
    [val,idx]=max(T);
    %find indeces of FITLEVELS
    k=idx;
    while T(k)>max(FITLEVELS);k=k+1;end;h=k-1;
    k=h;
    while T(k)>min(FITLEVELS);k=k+1;end;l=k+1;
    %exponential fit of transient decay
    XDATA=[h:8000]';
    FITX=[1:size(XDATA)]';
    FITY=T(XDATA);
    STARTP=[max(FITLEVELS),1000];
    OPTIONS=optimset('display','final','MaxFunEvals',2e6,'MaxIter',2e3,'TolFun',1e-9,'TolX',1e-8);
    P=lsqcurvefit('uptakefitfunction',STARTP,FITX,FITY,[STARTP(1)-0.02 100],[STARTP(1)+0.02 inf],OPTIONS);
    FT=uptakefitfunction(P,FITX);
    %plot result
    figure
    plot(FITX,FITY,'k');hold on
    plot(FITX,FT,'r')
    %find index of threshold decay
    k=1;
    while FT(k)>threshold && k<FITX(end);k=k+1;end;tindex=k-1;
    plot(tindex,FT(tindex),'b+');
    INDEX95(i)=tindex+XDATA(1)-1;
end
%% prepare data for fit and calculate initial conditions
timepoints=6000;
meanidx=round(mean(INDEX95));
TIME=[1:timepoints]'/scanrate;%time basis (s)
Z=zeros(size(TIME,1),size(STRANS,2));
Z0=zeros(size(STRANS,2),1);
V0=zeros(size(STRANS,2),1);
for i=1:size(STRANS,2);
    %relative coordinates z=(x-a)/a
    rest=mean(STRANS(end-500,i),1);
    Z(:,i)=(STRANS(meanidx:meanidx+timepoints-1,i)-rest)/rest;
    Z0(i)=Z(1,i);%initial relative sarcomere length
    %initial velocity
    %linear fit to the first 500 points in Z
    linfitpoints=100;
    XDATA=TIME(1:linfitpoints);
    YDATA=Z(1:linfitpoints,i);
    pz=polyfit(XDATA,YDATA,1);
    FITYDATA=polyval(pz,XDATA);
    V0(i)=scanrate*(FITYDATA(end)-FITYDATA(1))/linfitpoints;%initial velocity (1/s)
    %display fit result
    figure
    plot(TIME(1:500),Z(1:500,i),'k');hold on
    plot(XDATA,YDATA,'r');
    plot(XDATA,FITYDATA,'m');
end
Z(:,end-1:end)=[];
Z0(end-1:end)=[];
V0(end-1:end)=[];
%% fit sarcomere data
clear Z1 V1 FITZ FITP;global Z1 V1
STARTP=[-3,50,-0.01,-0.005];
%i=1;
for i=1:size(Z,2)
Z1=Z0(1);V1=V0(1);
%exp(-(X-P(3))*(P(1)-V1)/Z1).*(P(1)/P(2)*sin(P(2)*(X-P(3)))+Z1*cos(P(2)*(X-P(3))))+P(4);
OPTIONS=optimset('display','final','MaxFunEvals',2e6,'MaxIter',2e3,'TolFun',1e-5,'TolX',1e-5);
P=lsqcurvefit('relaxfitfunction',STARTP,TIME,Z(:,i),[-inf 10 -0.05 -0.001],[inf 100 0.1 0.001],OPTIONS)
FITZ(:,i)=relaxfitfunction(P,TIME);
FITP(i,:)=P;
%residuals
%plot fit result
figure
plot(TIME,Z(:,i));hold on
plot(TIME,FITZ(:,i),'r');
end
%% recover coefficients
clear B K
for i=1:size(FITP,1)
    P1=FITP(i,1);P2=FITP(i,2);Z1=Z0(i);V1=V0(i);
    B(i)=2/Z1*(P1-V1);
    K(i)=REST(i)*1.0e-9*(P2^2+((P1-V1)/Z1)^2);
end
    