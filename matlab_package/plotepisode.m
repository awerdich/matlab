%% prepare voltage traces
if exist('NORMDATA','var')==1
    FITPIXELS=zeros(size(NORMDATA,3),size(PIXELCOORDINATES,1));
    if exist('FILTERDATA','var')==1, FITPIXELSB=zeros(size(FILTERDATA,3),size(PIXELCOORDINATES,1)); end
    for i=1:size(PIXELCOORDINATES,1)
        FITPIXELS(:,i)=squeeze(NORMDATA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
        if exist('FILTERDATA','var')==1 && isempty(FILTERDATA)==0
            FITPIXELSB(:,i)=squeeze(FILTERDATA(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));   
            FITPIXELSB(:,i)=FITPIXELSB(:,i)-mean(FITPIXELSB(:,i));
            FITPIXELSB(:,i)=-FITPIXELSB(:,i)/max(FITPIXELSB(:,i));
        end
    end
end

%% prepare ca traces
if exist('EPISODER','var')
    FITPIXELS=zeros(size(EPISODER,3),size(PIXELCOORDINATES,1));
    for i=1:size(PIXELCOORDINATES,1)
        R=squeeze(EPISODER(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
        B=squeeze(BSCANS(PIXELCOORDINATES(i,2),PIXELCOORDINATES(i,1),:));
        %normalize trace
        R=R-mean(R(B));       
        FITPIXELS(:,i)=R/max(R);
    end
end

%% plot traces
figure
plot(FITPIXELS,'k')
hold on
plot(UPSTROKERANGE,FITPIXELS(UPSTROKERANGE,:),'b');
if exist('EPISODESTIM','var')
    plot(EPISODESTIM,'r')
end
hold off
%plot unnormalized data if available
if exist('FITPIXELSB','var')==1
    figure
    plot(FITPIXELSB);
end
    
    