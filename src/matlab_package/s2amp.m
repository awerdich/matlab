%this function measures the relative action potential amplitudes after an S2 stimulus
%pulse that follows an S1 pulse
%UPSTROKERANGE selected in normalizestack script must include S1 pulse
function [S2AMPMATRIX,S2AMPMATRIXBIN]=s2amp(NORMDATA,NSTIM,UPSTROKERANGE,scanrate)
%% define S1 and S2 timings
S2AMPMATRIX=zeros(size(NORMDATA,1),size(NORMDATA,2));
threshold=0.7;%stimulus detection threshold
STIMPEAK=[];
i=UPSTROKERANGE(1);
while i<length(NSTIM)
    if threshold<NSTIM(i)
        STIMPEAK=[STIMPEAK;i];
        while threshold<NSTIM(i) && i<length(NSTIM);i=i+1;end
    end
    i=i+1;
end
%% find S1 and S2 pulses and plot results
figure
plot(NSTIM,'k')
hold on
S1=[];S2=[];
for i=1:2:length(STIMPEAK)
    S1=[S1;STIMPEAK(i)];
    plot(STIMPEAK(i),NSTIM(STIMPEAK(i)),'r+')
    fprintf(['S1 at:',num2str(S1(end)),' (RED)\n']);
    if i+1<=length(STIMPEAK)
        S2=[S2;STIMPEAK(i+1)];
        plot(STIMPEAK(i+1),NSTIM(STIMPEAK(i+1)),'b+')
        fprintf(['S2 at:',num2str(S2(end)),' (BLUE)\n']);
    end
end
%ask for stimulus rate if only one S1 pulse in recording
if length(S1)<2
    s1interval=input('S1 interval (ms):');
    S1(2,1)=S1(1)+s1interval*1.0e-3*scanrate;
end    
%print coupling interval
fprintf(['S1-S2 coupling interval (ms):',num2str((S2(1)-S1(1))/scanrate*1.0e3),'\n']);
%% measure amplitudes for all nonzero pixels
minrisetime=50;%signals shorter than that are defined as artifacts
apthreshold=0.6;
bin=0.2;%divide the amplitudes in equal intervals
BININTERVAL=[];b=0;while b<1;BININTERVAL=[BININTERVAL;b];b=b+bin;end;BININTERVAL=[BININTERVAL;1];
S2AMPMATRIX=zeros(size(NORMDATA,1),size(NORMDATA,2));
S2AMPMATRIXBIN=zeros(size(S2AMPMATRIX));
ARTIFACTMATRIX=zeros(size(S2AMPMATRIX));
hdl = waitbar(0,['Measuring amplitudes. Please wait.']);
for i=1:size(NORMDATA,1)
    for j=1:size(NORMDATA,2)
        PIXEL=squeeze(NORMDATA(i,j,:));
        if abs(mean(PIXEL(:)))>0
            
            %Measure amplitudes of membrane response to S1 and S2 pulses
            S1AMP=[];%AP amplitudes following S1 pulses go in here
            S2AMP=[];%AP amplitudes sollowing S2 pulses go in here
            %measure amplitudes for S1 pulses
            for k=1:length(S1)
                %k is the scan of an S1 pulse
                %find the upstroke and the downstroke phases
                upscan=0;downscan=0;
                n=S1(k);while PIXEL(n)<apthreshold && n<length(PIXEL);n=n+1;end;upscan=n+1;
                %go a little higher to make sure that the next scan is not
                %below apthreshold
                while PIXEL(n)<apthreshold+0.1 && n<length(PIXEL);n=n+1;end;
                while apthreshold<PIXEL(n) && n<length(PIXEL);n=n+1;end;downscan=n-1;
                %try to find the maximum of the PIXEL time series between
                %upscan and downscan
                if upscan>0 && downscan>0
                    [maxval,maxidx]=max(PIXEL(upscan:downscan));maxidx=maxidx+upscan-1;
                    S1AMP=[S1AMP;mean(PIXEL(maxidx-10:maxidx+10))];
                end
            end
            %measure amplitude for S2 pulse following the first S1 pulse
            %it should be between S2 and the second S1
            AMP=[];
            for s=S2(1):S1(2)
                AMP=[AMP;[s,PIXEL(s)]];    
            end
            %sort amplitudes
            SAMP=sortrows(AMP,2);
            %average the top 10% of amplitudes
            avscans=round(0.1*length(SAMP));
            s2amp=mean(SAMP(end-avscans:end,2))/mean(S1AMP);
            S2AMPMATRIX(i,j)=s2amp;
                  
            %determine binned S2AMPMATRIX
            %if amp<second bin, count it as 0
            if s2amp<=BININTERVAL(2)
                S2AMPMATRIXBIN(i,j)=0.1;
            %if amp>last bin, count it as 1
            elseif BININTERVAL(end-1)<s2amp
                S2AMPMATRIXBIN(i,j)=1;
            end
            %if 20%<s2amp<80% then sort them into bin-size intervals
            for b=2:length(BININTERVAL)-2
                if BININTERVAL(b)<s2amp && s2amp<=BININTERVAL(b)+bin/2
                    S2AMPMATRIXBIN(i,j)=BININTERVAL(b);
                elseif BININTERVAL(b)+bin/2<s2amp && s2amp<=BININTERVAL(b+1)
                    S2AMPMATRIXBIN(i,j)=BININTERVAL(b+1);
                end
            end
            
            %Find artifacts by measuring amplitude and width of response
            B=[1:avscans];%baseline
            OAMP=AMP(:,2)-mean(SAMP(B,2));
            NAMP=OAMP/max(OAMP(:));
            %find rise time of NAMP
            [maxval,maxidx]=max(NAMP);
            left=maxidx;while 0.5<NAMP(left) && 1<left;left=left-1;end
            right=maxidx;while 0.5<NAMP(right) && right<length(NAMP);right=right+1;end
            %criteria for artifact: amplitude>20% and rise time<20 ms
            if BININTERVAL(2)<s2amp && ((left==1 && right<length(NAMP)) || (right-left)<round(minrisetime*1.0e-3*scanrate))
                ARTIFACTMATRIX(i,j)=1;
            end  
        end
    waitbar(i/size(NORMDATA,1));   
    end
end
close(hdl)
%apply artifact detection to binned matrix
 S2AMPMATRIXBIN(ARTIFACTMATRIX==1)=0.1;
end
                    
                    
                    
                    
                
                
                
        
        