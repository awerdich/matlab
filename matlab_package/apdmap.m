%% User parameter
DATA=NORMDATA;
SIGNALS=ASIGNALPIXELS;
apdthreshold=0.2;
if exist('SCATTER','var')==0;SCATTER=ones(size(ASIGNALPIXELS));end
SIGNALS(SCATTER<0.5)=0;
APDMATRIX=zeros(size(ASIGNALPIXELS));
VMAXMATRIX=zeros(size(ASIGNALPIXELS));%maximum upstroke velocities of NORMDATA
VMAXMATRIXSM=zeros(size(ASIGNALPIXELS));%maximum upstroke velocities of smoothed NORMDATA
golaywindow=31;%smoothing filter window sent to CAM on 10/6/2010 for cplexapd
tstart=tic;%start timer
hdl = waitbar(0,['calculating APDs and upstroke velocities']);
for i=1:size(DATA,1)
    parfor j=1:size(DATA,2)
    %% calculate APD for pixel (i,j)                     
        if ASIGNALPIXELS(i,j)==1
            PIXEL=squeeze(DATA(i,j,:));
            BASELINE=squeeze(BASELINESCANS(i,j,:));
            [upidx,downidx,T]=apd(PIXEL,apdthreshold,UPSTROKERANGE,BASELINE,scanrate);
            if upidx>0 && downidx>0
                APDMATRIX(i,j)=(downidx-upidx)*1000/scanrate;
                %save unfiltered upstroke velocity
                VMAXMATRIX(i,j)=max(diff3(T(UPSTROKERANGE)))*scanrate;
            end
            if upidx>0
                %apply smoothing algorithm to calculate first derivative
                [D0,D1,D2]=smoothdiff(T(UPSTROKERANGE),golaywindow);
                %find maximum of first differential
                VMAXMATRIXSM(i,j)=max(D1)*scanrate;
            end
        end
%% start over with new pixel 
    end
    waitbar(i/size(DATA,1));
end
close(hdl)
elapsedtime=toc(tstart);fprintf(['Elapsed time: ',num2str(elapsedtime),' s\n']);
%% save data
un='[APD]=ms;[Vmax]=1/s';
if exist('datapath')==0 || isempty(datapath)==1 || length(datapath)<4
    if exist('stackpath')==0
        path='D:\';
    end
    datapath=uigetdir(path,'DATA FOLDER'); 
    datapath=[datapath,'\'];
end

apdfile=[normfile(1:end-5),'APD',num2str(100-apdthreshold*100),'.mat'];

save([datapath,apdfile],'APDMATRIX','VMAXMATRIX','VMAXMATRIXSM','apdthreshold','golaywindow','BASELINE','UPSTROKERANGE','datapath','scanrate');