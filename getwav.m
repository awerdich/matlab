regionpol_v = 'v_POL';
wavefronts_v = calculatew(DATABASE, regionpol_v);
wav = wavefronts_v(:,3)


%%

function WAVEFRONTS = calculatew(DATABASE, regionpol)
%% LOOP
% find activation times for pixels that are far away
WAVEFRONTS=[];SET=[];
for id=1:length(DATABASE)
    %% extract activation map from database
    scanrate=DATABASE(id).scanrate;%scans/s
    if isfield(DATABASE,'pixelcalfactorj')==1
        pixelcalfactor=DATABASE(id).pixelcalfactorj;%[um/pixel]
    elseif isfield(DATABASE,'pixecalfactor')==1
        pixelcalfactor=DATABASE(id).pixelcalfactor;
    end
    if isempty(DATABASE(id).(regionpol))==0
        SET=[SET;id];%dataset used
    %activation matrix
    POL=DATABASE(id).(regionpol);
    %% normalize 
    %determine data range
    T=[];
    for i=1:size(POL,1)
        for j=1:size(POL,2)
            if POL(i,j)>0
                T=[T;POL(i,j)];
            end
        end
    end
    tmin=min(T);tmax=max(T);%activation time data limit
    lowlimit=0;highlimit=1.0;%additional limit for lowest and highest activation times
    
    %re-map activation times
    NPOL=zeros(size(POL));
    LNPOL=zeros(size(POL));
    for i=1:size(POL,1)
        for j=1:size(POL,2)
            if POL(i,j)>0
                    %map entire activation time into [0,1]
                    NPOL(i,j)=(POL(i,j)-tmin)/(tmax-tmin);
                    %map activation times into [lowlimit,highlimit]
                    if NPOL(i,j)<lowlimit
                        LNPOL(i,j)=0;
                    elseif lowlimit<=NPOL(i,j) && NPOL(i,j)<highlimit
                        LNPOL(i,j)=(NPOL(i,j)-lowlimit)/(highlimit-lowlimit);
                    else
                        LNPOL(i,j)=1;
                    end
            end
        end
    end
    
    %Re-size data matrix
    resizefactor=2;%resize divergence image for plotting
    %signal matrix
    SIGNAL=zeros(size(POL));
    SIGNAL(POL>0)=1;
    %re-size data
    RNPOL=imresize(LNPOL,resizefactor,'bicubic');
    RSIGNAL=imresize(SIGNAL,resizefactor,'bicubic');
    RSIGNAL(RSIGNAL<0.9)=0;RSIGNAL(RSIGNAL>=0.9)=1;
    %fix interpolated activation matrix at signalk boundary
    RNPOL(RSIGNAL==0)=0;
    
    %thresholding 
    dt=0.1;%time interval for thresholding
    t=dt;%start time
    n=floor(1/dt);%number of time steps
    TPOL=zeros(size(RNPOL,1),size(RNPOL,2),n);
    ISLANDS=[];%number of wavefronts
    while t+dt<=1
        TINTERVAL=[t;t+dt];
        TPOL=zeros(size(RNPOL));
        for i=1:size(RNPOL,1)
            for j=1:size(RNPOL,2)
                %exclude pixels that are off
                if RSIGNAL(i,j)==1
                    if TINTERVAL(1)<=RNPOL(i,j) && RNPOL(i,j)<=TINTERVAL(2)
                        TPOL(i,j)=1;
                    end
                end
            end
        end
        %calculate the connectivity of TPOL
        CC=bwconncomp(TPOL,8);
        ISLANDS=[ISLANDS;[id,t,CC.NumObjects]];%islands as 
        t=t+dt;
    end
    %maximum number of islands
    [val,idx]=max(ISLANDS(:,3));
    WAVEFRONTS=[WAVEFRONTS;ISLANDS(idx,:)];
    save WAVEFRONTS.txt WAVEFRONTS -ascii -tabs
    end
end
end