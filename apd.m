%function to calculate APD in scans
function [UPIDX,DOWNIDX,MEANTRANS]=apd(DATA,threshold,UPSTROKERANGE,BASELINE,scanrate)
UPIDX=zeros(size(DATA,2),1);
DOWNIDX=zeros(size(DATA,2),1);
if threshold<0.3
    UFITINTERVAL=[threshold-0.1;threshold+0.6];%fit interval in upstroke
    DFITINTERVAL=[threshold-0.1;threshold+0.4];%fit interval in decay
else
    UFITINTERVAL=[threshold-0.3;threshold+0.3];%fit interval in upstroke
    DFITINTERVAL=[threshold-0.3;threshold+0.3];%fit interval in decay
end
NORMTRANS=zeros(size(DATA));%collect normalized transients
maxscans=round(5*1e-3*scanrate)+1;%number of scans to average for measuring the maxiumum;

for i=1:size(DATA,2) 
            PIXEL=DATA(:,i);
            %correct baseline
            OPIXEL=PIXEL-mean(PIXEL(BASELINE));
            %define upstroke range as the scan interval PIXEL([min,max])
            %with min<max
            RANGE=UPSTROKERANGE;
            [minval,minidx]=min(OPIXEL(RANGE));minidx=minidx+RANGE(1)-1;
            [maxval,maxidx]=max(OPIXEL(RANGE));maxidx=maxidx+RANGE(1)-1;
            %shorten RANGE until minidx<maxidx
            while maxidx<minidx && length(RANGE)>1
                RANGE=RANGE(1:end-1);%
                [minval,minidx]=min(OPIXEL(RANGE));minidx=minidx+RANGE(1)-1;
                [maxval,maxidx]=max(OPIXEL(RANGE));maxidx=maxidx+RANGE(1)-1;
            end
            %upstroke range is [minidx,maxidx]          
            maximum=mean(OPIXEL(maxidx-ceil(maxscans/2):maxidx+ceil(maxscans/2)));
            NPIXEL=OPIXEL/maximum;
            NORMTRANS(:,i)=NPIXEL;
            
            %search for threshold level in the upstroke starting from max
            up=maxidx;while NPIXEL(up)>threshold && up>1;up=up-1;end
            %up is now within UFITINTERVAL per definition
            
            if up>RANGE(1)
                %find the scans in the UFITINTERVAL by starting at the
                %maximum and going backwards
                upfitlow=up;while NPIXEL(upfitlow)>UFITINTERVAL(1) && upfitlow>1;upfitlow=upfitlow-1;end
                if upfitlow>1
                    upfithigh=up;while NPIXEL(upfithigh)<UFITINTERVAL(2);upfithigh=upfithigh+1;end
                    UPXDATA=[upfitlow:upfithigh]';UPYDATA=NPIXEL(UPXDATA);
            
                    %fit
                    if length(UPXDATA)>10
                        fitorder=5;
                    else
                        fitorder=1;
                    end
            
                    [upfit,S,MU]=polyfit(UPXDATA,UPYDATA,fitorder);
                    %calculate threshold index from fitted polynom in finer
                    UPXDATAFIT=[UPXDATA(1):0.001:UPXDATA(end)]';
                    UF=(UPXDATAFIT-MU(1))./MU(2);%coordinate transformation to calculate fitted polynom
                    UPYDATAFIT=polyval(upfit,UF);
                    %search for new threshold on fit
                    upfitidx=1;while UPYDATAFIT(upfitidx)<threshold && upfitidx<length(UPYDATAFIT);upfitidx=upfitidx+1;end
                    UPIDX(i)=UPXDATAFIT(upfitidx);
                else
                    UPIDX(i)=0;
                end
            else
                UPIDX(i)=0;
            end
            
            %repolarization phase
            %start searching from maximum
            down=maxidx;while NPIXEL(down)>threshold && down<length(NPIXEL);down=down+1;end
            %down is now within DFITINTERVAL per definition
            if UPIDX(i)>0 && down<(length(NPIXEL)) %-length(RANGE))
                %continue to search downstroke only if upstroke was found
                %and threshold level found within UPSTROKERANGE
                %find fit interval boundaries: start with high border
                downfithigh=down;while NPIXEL(downfithigh)<DFITINTERVAL(2) && downfithigh>maxidx;downfithigh=downfithigh-1;end
                downfitlow=down;while NPIXEL(downfitlow)>DFITINTERVAL(1) && downfitlow<length(NPIXEL);downfitlow=downfitlow+1;end
                DOWNXDATA=[downfithigh:downfitlow]';DOWNYDATA=NPIXEL(DOWNXDATA);
                %fit
            
                if length(DOWNXDATA)>10
                    fitorder=5;
                else
                    fitorder=1;
                end
                            
                [downfit,S,MU]=polyfit(DOWNXDATA,DOWNYDATA,fitorder);
                %calculate threshold index from fitted polynom in finer
                DOWNXDATAFIT=[DOWNXDATA(1):0.001:DOWNXDATA(end)]';
                DF=(DOWNXDATAFIT-MU(1))./MU(2);%coordinate transformation to calculate fitted polynom
                DOWNYDATAFIT=polyval(downfit,DF);
                %search for new index on fit
                downfitidx=1;while DOWNYDATAFIT(downfitidx)>threshold && downfitidx<length(DOWNYDATAFIT);downfitidx=downfitidx+1;end
                DOWNIDX(i)=DOWNXDATAFIT(downfitidx);
            else
                DOWNIDX(i)=0;
            end
end
MEANTRANS=mean(NORMTRANS,2);%return normalized mean transient
        