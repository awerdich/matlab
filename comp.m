%% get original ratios and calibration constants
ntrace=(FrameYLim(end)-FrameYLim(1)+1)*(FrameXLim(end)-FrameXLim(1)+1);
TRACE=[];MINR=[];MAXR=[];K=[];F=[];
for i=FrameYLim(1):FrameYLim(end)
    for j=FrameXLim(1):FrameXLim(end)
        if SIGNALS(i,j)>0
            TRACE=[TRACE,squeeze(TRACEDATA(i,j,:))];
            MINR=[MINR;RMIN(i,j)];
            MAXR=[MAXR;RMAX(i,j)];
            K=[K;KD(i,j)];
            F=[F;F0(i,j)];
        end
    end
end
if exist('EPISODESTIM')==0
    EPISODESTIM=[];
end
R=[mean(TRACE,2)];
rmin=mean(MINR);
rmax=mean(MAXR);
k=mean(K);
f=mean(F);
%% increase R to 90% of Rmax
%max(R)/Rmax=fract
FRACT=[0.2:0.2:0.8];
BASELINE=[239:255];offset=mean(R(BASELINE));
%offset
R0=R-offset;    
CA1=zeros(size(R,1),length(FRACT));
NCA1=zeros(size(CA1));

%normal Ca trace
CA=k*f*(R-rmin)./(rmax-R);
NCA=CA-mean(CA(BASELINE));
NCA=NCA/max(NCA(1:200));

for j=1:length(FRACT)
    for i=1:length(TRACE)
        %re-normalize R
        newmax=FRACT(j)*rmax-offset;
        R1=R0/max(R0)*newmax;
        R1=R1+offset;
        %calibrate CA1
        CA1(i,j)=k*f*(R1(i)-rmin)/(rmax-R1(i));
    end
    %normalize
    NCA1(:,j)=CA1(:,j)-mean(CA1(BASELINE,j));
    NCA1(:,j)=NCA1(:,j)/max(NCA1(1:200,j));
end