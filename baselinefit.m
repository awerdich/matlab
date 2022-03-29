function FITBASELINE=baselinefit(PIXEL,XDATA,METHOD)
%fits baseline of PIXEL trace to defined model

switch METHOD
    case 'linear'
        %polynomial fit model
        FITX=XDATA;FITY=PIXEL(XDATA);
        P=polyfit(XDATA,PIXEL(XDATA),1);
        FITBASELINE=polyval(P,1:length(PIXEL))'; 
    case 'exp'
        %exponential fit model
        %sort XDATA
        SXDATA=sortrows(XDATA,1);
        %remove offset
        offset=mean(PIXEL(SXDATA(1:20)));
        OPIXEL=PIXEL-offset;
        %normalize 
        if abs(mean(PIXEL))>0
            maxdata=max(OPIXEL(SXDATA));
            NPIXEL=OPIXEL/maxdata;
        else
            maxdata=0;
            NPIXEL=PIXEL;
        end
        
        %define fit data
        FITX=SXDATA;
        FITY=NPIXEL(SXDATA);
    
        %startparameter
        startx=length(PIXEL)/2;
        starta=-startx/log(0.5);
    
        foptions=fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,-inf],...
            'Upper',[inf,inf],...
            'StartPoint',[starta,0],...
            'TolFun',1.0e-3,...
            'MaxIter',50,...
            'MaxFunEvals',50);
    
        fitfunc=fittype('1-exp(x0-x/a)','options',foptions);
    
        [y,p]=fit(FITX,FITY,fitfunc);
        FITDATA=feval(y,[1:length(PIXEL)]');
        FITBASELINE=(FITDATA*maxdata)+offset;     
end
 