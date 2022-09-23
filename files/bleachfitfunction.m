function FITDATA=bleachfitfunction(DATA)
%% define baseline 
T=DATA(:,1);
% plot trace T
figure
plot(T,'k');
%select range
fprintf('select fit range\n')
[X,Y]=ginput(2);
RANGE=[round(min(X)):round(max(X))];
TR=T(RANGE);
mainfig=figure;
plot(TR,'k')
%enter baseline
fprintf('define baseline intervals. ESC to end\n');
XPOINTS=[];
button=1;
while button~=27
    [xi,yi,button]=ginput(1);
    if button~=27
        XPOINTS=[XPOINTS;round(xi)];
    end
end
%mark baseline in mainfig
XDATA=[];INTERVALX=[];
figure(mainfig);hold on
for i=1:2:length(XPOINTS)-1
    INTERVALX=(XPOINTS(i):1:XPOINTS(i+1))';  
    XDATA=[XDATA;INTERVALX];
    plot(INTERVALX,TR(INTERVALX,:),'m')
end

%% fit all traces
XDATA=sort(XDATA);
FITDATA=zeros(length(RANGE),size(DATA,2));
for i=1:size(DATA,2)
    T=DATA(RANGE,i);
    P=fit(XDATA,T(XDATA),'poly1');
    %correct data
    FT=T-P(RANGE);
    %correct offset
    OFT=FT-mean(FT(XDATA));
    FITDATA(:,i)=OFT;
end



