clear
%% load data file
if exist('loadpath','var')==0 || length(loadpath)<2
    [datafile, loadpath] = uigetfile('*.mat', 'Pick a DATABASE file');
end
load([loadpath,datafile]);
%% field definitions
fprintf(['DATABASE:',num2str(datafile),'\n']);
region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
prop=input('impulse initiation 0=atrium 1=av 2=ventricle 3=fstim 4=any:');

if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
|| strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
|| strcmp(region,'vic')==1 || strcmp(region,'voc')==1

   %CV FIELD DEFINITIONS 
   %regionangleunitvectorxy=[region,'_ANGLEUNITVECTORXY'];%angle unit vector defined for this region
   regionmaxrmse=[region,'_maxrmse'];%maxRMSE
   regionvelinterval=[region,'_VELINTERVAL'];%velocity inteval for fitted vectors
   regionv=[region,'_vel'];%mean propagation velocity
   regionvstd=[region,'_velstd'];%standard deviation of FINAL vector list (SORTV)
   regionsignalframe=[region,'_SIGNALFRAME'];
   regionframearea=[region,'_framearea_um'];%active tissue area in frame (um^2, from SIGNALFRAME)
   regiontotalarea=[region,'_totalarea_um'];%estimated total tissue area (um^2, from SIGNAL)
   regionfits=[region,'_numfits'];%total number fits performed in frame
   regionwindowx=[region,'_windowx'];%frameumwidth
   regionwindowy=[region,'_windowy'];%frameumheight
   regionwindowcenter=[region,'_FRAMECENTER'];%frame center location
   regionframevi=[region,'_VFRAMEi'];%conduction velocity vector i-component [pixel/frame]
   regionframevj=[region,'_VFRAMEj'];%conduction velocity vector j-component [pixel/frame]
   regionframex=[region,'_XFRAMEINTERP0'];%x-coordinates of frame region
   regionframey=[region,'_YFRAMEINTERP0'];%y-coordinates of frame region
   regionsignalframeinterp=[region,'_SIGNALFRAMEINTERP0'];%Binary image: Signals in frame
   regionplotbin=[region,'_PLOTVELOCITYBIN'];%VPLOTBIN binary image with filtered signals
   regionsignalframeinterpolygon=[region,'_POLYROI'];%Binary image: Signals in polygon
   regionrmse=[region,'_RMSE'];%RMSE for each vector in frame
   regionpol=[region,'_POL'];
   regionrepol=[region,'_REPOL'];

   %APD FIELD DEFINITIONS
   regionframedata=[region,'_TRACEDATA'];%ALL TRANSIENTS IN FRAME (i,j,:)
   regionmeantrace=[region,'_TRACE'];%MEAN FRAME TRANSIENT
   regionapdlist=[region,'_APDLIST_ms'];%list of APDs in region (ms)
   apdthreshold=0.2;
   regionapd=[region,'_meanapd_ms',num2str(apdthreshold*100)];%AVERAGE TRANSIENT DURATION (ms) 
   regionvmaxlist=[region,'_VMAXLIST_s'];%list of VMAX values in region (1/s);
   regionvmaxlistsm=[region,'_VMAXLISTSM_s'];%list of VMAX values from smoothed transients (1/s);
   regionvmax=[region,'_meanvmax_s'];%average maximum upstroke velocity (1/s);
   regionvmaxsm=[region,'_meanvmaxsm_s'];%average vmax from smoothed transients(1/s)
   regionstimtrace=[region,'_STIMTRACE'];%stimulus transient;
else
    fprintf('Wrong input! \n');
    return
end

%COLLECT DATA
%smmary definitions

%number of entry in database
ENTRY=[];ENTRYAPD=[];

%rates
RATE=[];RATEAPD=[];

%velocities
VEL=[];VELSTD=[];

%Activation times
POLROI=[];REPOLROI=[];
POLROIN=[];

%divergence
DIV=[];
DIVOUTLIERS=[];
DIVPOSITIVE=[];
DIVNEGATIVE=[];

%angles
ANG=[];ANGSTD=[];ANGLISTN=[];

%APD
APD=[];APDSTD=[];APDSEM=[];

%ap upstroke
VMAX=[];VMAXSTD=[];VMAXSEM=[];
VMAXSM=[];VMAXSTDSM=[];VMAXSEMSM=[];

%data collection
 hdl = waitbar(0,'PROCESSING DATABASE ENTRIES');
for id=1:length(DATABASE)
    %collect conduction data
    if (prop==DATABASE(id).sinusloc || prop==4)&& isfield(DATABASE,(regionv))...
            && isempty(DATABASE(id).(regionv))==0   
        %heartrates
        ENTRY=[ENTRY;id];
        RATE=[RATE;DATABASE(id).rate];
        
       %VELOCITY VECTORS 
        VPLOTBIN=DATABASE(id).(regionplotbin);%Binary image with filtered velocities
        if isfield(DATABASE,(regionsignalframeinterpolygon))==1 && isempty(DATABASE(id).(regionsignalframeinterpolygon))==0
            POLYROI=DATABASE(id).(regionsignalframeinterpolygon);%Binar image with ROIs
            fprintf(['processing polygon ROI for dataset:',num2str(id),'/',num2str(length(DATABASE)),'\n']);
        else
            clear POLYROI;%remove POLYROI if it existed in the previous dataset
            fprintf(['processing rectangular ROI for dataset:',num2str(id),'/',num2str(length(DATABASE)),'\n']);
        end
        RMSE=DATABASE(id).(regionrmse);
        maxrmse=DATABASE(id).(regionmaxrmse);%marmse
        scanrate=DATABASE(id).scanrate;%[frames/s]
        if isfield(DATABASE,'pixelcalfactori')==1 && isempty(DATABASE(id).pixelcalfactori)==0
            pixelcalfactor_i=DATABASE(id).pixelcalfactori;%[um/pixel]
            pixelcalfactor_j=DATABASE(id).pixelcalfactorj;%[um/pixel]
        else
            pixelcalfactor_i=DATABASE(id).pixelcalfactor;%[um/pixel]
            pixelcalfactor_j=DATABASE(id).pixelcalfactor;%[um/pixel]
        end

        %% calculate divergence for current ROI
        %velocity components for ROI
        MI=DATABASE(id).(regionframevi);%[pixel/frame]
        MJ=DATABASE(id).(regionframevj);%[pixel/frame]
        %divergence for ROI
        [DIVFIELD,VOUTLIERS]=cvdivergence(MJ,MI);
        %collect velocities and divergences within ROI
        VECTORLIST=[];
        DIVLIST=[];%all divergence values
        DIVLISTPOSITIVE=[];%only positive divergence values
        DIVLISTNEGATIVE=[];%only negative divergence values
        for i=1:size(VPLOTBIN,1)
            for j=1:size(VPLOTBIN,2)
                %exclude filtered pixels
                if VPLOTBIN(i,j)==1 && (0<RMSE(i,j) && RMSE(i,j)<maxrmse)
                    %absolute velocity at (i,j)
                    vi=MI(i,j)*scanrate*pixelcalfactor_i*1.0e-3;%[pixel/frame*frame/s*um/pixel*1e-3mm/um]
                    vj=MJ(i,j)*scanrate*pixelcalfactor_j*1.0e-3;%[mm/s]
                    if exist('POLYROI','var')==1
                        %collect velocities within POLYROI
                        if POLYROI(i,j)==1    
                            VECTORLIST=[VECTORLIST;[i,j,vi,vj,norm([vi,vj],2)]];
                            %exclude bad pixels (NaN) and pixels with velocity outliers
                            if isnan(DIVFIELD(i,j))==0 && DIVFIELD(i,j)~=0 && VOUTLIERS(i,j)==0  
                                DIVLIST=[DIVLIST;[i,j,DIVFIELD(i,j)*scanrate]];%[1/s]
                                if DIVFIELD(i,j)>0
                                    DIVLISTPOSITIVE=[DIVLISTPOSITIVE;[i,j,DIVFIELD(i,j)*scanrate]];
                                else
                                    DIVLISTNEGATIVE=[DIVLISTNEGATIVE;[i,j,DIVFIELD(i,j)*scanrate]];
                                end
                            end
                        end
                    else
                        %add vectors in rectangulat ROI without POLYNOM
                        VECTORLIST=[VECTORLIST;[i,j,vi,vj,norm([vi,vj],2)]];
                        if isnan(DIVFIELD(i,j))==0 && DIVFIELD(i,j)~=0 && VOUTLIERS(i,j)==0  
                                DIVLIST=[DIVLIST;[i,j,DIVFIELD(i,j)*scanrate]];%[1/s]
                                if DIVFIELD(i,j)>0
                                    DIVLISTPOSITIVE=[DIVLISTPOSITIVE;[i,j,DIVFIELD(i,j)*scanrate]];
                                else
                                    DIVLISTNEGATIVE=[DIVLISTNEGATIVE;[i,j,DIVFIELD(i,j)*scanrate]];
                                end
                        end
                    end
                end
            end
        end
        
        
        %% calculate angles
        %ANGLES
        %calculate mean velocity vector and angles
        MEANVECTORij=[sum(VECTORLIST(:,3));sum(VECTORLIST(:,4))]/size(VECTORLIST,1);
        if isfield(DATABASE,'ANGLEUNITVECTORXY') && isempty(DATABASE(id).ANGLEUNITVECTORXY)==0
            REFERENCEVECTORxy=DATABASE(id).ANGLEUNITVECTORXY;
        else
            REFERENCEVECTORxy=DATABASE(id).angleunitvectorxy;
        end
        REFERENCEVECTORij=[REFERENCEVECTORxy(2);REFERENCEVECTORxy(1)];
        ANGLIST=zeros(size(VECTORLIST,1),2);
        for i=1:size(VECTORLIST,1)
            VELOCITYVECTORij=[VECTORLIST(i,3);VECTORLIST(i,4)];
            %ANGLIST(i,:)=[vectorangle(MEANVECTORij,VELOCITYVECTORij),VECTORLIST(i,5)];
            ANGLIST(i,:)=[vectorangle(REFERENCEVECTORij,VELOCITYVECTORij),VECTORLIST(i,5)];
        end
        %% collect results for current heart
        VEL=[VEL;mean(VECTORLIST(:,5))];
        VELSTD=[VELSTD;std(VECTORLIST(:,5))];
        ANGLISTN=[ANGLISTN;ANGLIST];
        DIV=[DIV;mean(DIVLIST(:,3))];
        DIVOUTLIERS=[DIVOUTLIERS;bwarea(VOUTLIERS)];%number of pixels excluded from divergence calculation for this id.
        DIVPOSITIVE=[DIVPOSITIVE;[mean(DIVLISTPOSITIVE(:,3))]];
        DIVNEGATIVE=[DIVNEGATIVE;[mean(DIVLISTNEGATIVE(:,3))]];
        
%         %calculate times to polarize entire ROI (IGNORE POLYNOM!)
%         PLIST=[];POL_ROI=[];SORTPLIST=[];
%         if isfield(DATABASE,regionpol)==1 && isempty(DATABASE(id).(regionpol))==0  
%             POL_ROI=DATABASE(id).(regionpol);
%             %FPOL=fdpolmatrix(POL_ROI,0,0);
%             FPOL=POL_ROI;
%             for i=1:size(FPOL,1)
%                 for j=1:size(FPOL,2)
%                     if FPOL(i,j)>0
%                         PLIST=[PLIST;FPOL(i,j)];
%                     end
%                 end
%             end
%             %sort polarization times
%             SORTPLIST=sort(PLIST);
%             SORTPLIST=SORTPLIST-min(PLIST);
%             %calculate the time difference
%             maxtime=mean(SORTPLIST(end-20:end-10));
%             DELTAPOLms=maxtime/scanrate*1.0e3;
%             POLROI=[POLROI;DELTAPOLms];
%         end
%         
        
        %calculate times to repolarize ROI
%         RLIST=[];
%         if isfield(DATABASE,regionrepol)==1
%             REPOL_ROI=DATABASE(id).(regionrepol);
%             %FREPOL=fdpolmatrix(REPOL_ROI,0,0);
%             FREPOL=REPOL_ROI;
%             for i=1:size(FREPOL,1)
%                 for j=1:size(FREPOL,2)
%                     if FREPOL(i,j)>0
%                         RLIST=[RLIST;FREPOL(i,j)];
%                     end
%                 end
%             end
%             %sort repolarization times
%             SORTRLIST=sort(RLIST);
%             SORTRLIST=SORTRLIST-min(RLIST);
%             maxtime=mean(SORTPLIST(end-20:end-10));
%             %calculate the time difference
%             DELTAREPOLms=maxtime/scanrate*1.0e3;
%             REPOLROI=[REPOLROI;DELTAREPOLms];
%         end
        
        %mean angle and standard deviation should be based on [0,180]
        ANGLIST180=ANGLIST;
        for i=1:size(ANGLIST,1)
            if ANGLIST(i,1)>180
                ANGLIST180(i,1)=360-ANGLIST(i,1);
            end
        end
        ANG=[ANG;mean(ANGLIST180(:,1))];
        ANGSTD=[ANGSTD;std(ANGLIST180(:,1))];
    end
    
    %APD data
    APDLIST=[];VMAXLIST=[];VMAXLISTSM=[];
    if prop==DATABASE(id).sinusloc || prop==4
    if DATABASE(id).rate>=0 && isfield(DATABASE,(regionapd)) && ...
        isempty(DATABASE(id).(regionapd))==0
    
    %entry
    ENTRYAPD=[ENTRYAPD;id];
    
    %rates
    RATEAPD=[RATEAPD;DATABASE(id).rate];
     
    %APD
    %extract APDs and TRACES in POLYROI (if defined)
    APDLISTFRAME=DATABASE(id).(regionapdlist);%APDLIST in FRAME with SIGNAL coordinates [i,j,APD(i,j)]
    TRACEDATAFRAME=DATABASE(id).(regionframedata);%APD traces in FRAME
    if isfield(DATABASE,(regionvmaxlistsm))==1 && isempty(DATABASE(id).(regionvmaxlistsm))==0
        VMAXLISTSMFRAME=DATABASE(id).(regionvmaxlistsm);
        VMAXLISTSM=[];
    end   
    if exist('POLYROI','var')==1
        %figure out which line is in APDLISTFRAME is activated POLYROI
        POLYINDEX=[];%line numbers that are within POLYROI
        XFRAMEINTERP_0=DATABASE(id).(regionframex);%absolute x coordinates of the FRAME
        YFRAMEINTERP_0=DATABASE(id).(regionframey);
        %check create new POLYROI based on APDs found
        APDFRAMEBIN=zeros(size(POLYROI));
        APDLIST=[];
        TRACEDATA=[];
        for i=1:length(APDLISTFRAME)
            x=APDLISTFRAME(i,2);%x-coordinate of APD line i
            y=APDLISTFRAME(i,1);%y-coordinate of APD line i
            %search for corresponding index in the X,Y FRAME
            jframe=1;while XFRAMEINTERP_0(1,jframe)~=x && jframe<size(XFRAMEINTERP_0,2);jframe=jframe+1;end
            iframe=1;while YFRAMEINTERP_0(iframe,1)~=y && iframe<size(YFRAMEINTERP_0,1);iframe=iframe+1;end
            %check if (iframe,jframe) are activated in POLYROI
            if POLYROI(iframe,jframe)==1
                APDLIST=[APDLIST;APDLISTFRAME(i,:)];
                TRACEDATA=[TRACEDATA,TRACEDATAFRAME(:,i)];
                APDFRAMEBIN(iframe,jframe)=1;
                if isfield(DATABASE,(regionvmaxlistsm))==1 && isempty(DATABASE(id).(regionvmaxlistsm))==0
                    VMAXLISTSM=[VMAXLISTSM;VMAXLISTSMFRAME(i,:)];
                end   
            end    
        end
    else
        APDLIST=APDLISTFRAME;
        TRACEDATA=TRACEDATAFRAME;
        if isfield(DATABASE,(regionvmaxlistsm))==1 && isempty(DATABASE(id).(regionvmaxlistsm))==0
            VMAXLISTSM=VMAXLISTSMFRAME;
        end   
    end
    
    APD=[APD;mean(APDLIST(:,3))];
    APDSTD=[APDSTD;std(APDLIST(:,3))];
    APDSEM=[APDSEM;std(APDLIST(:,3))/sqrt(size(APDLIST,1))];

    %upstroke velocities
    %calculate upstroke velocities from tracedata
    tracedatascanrate=DATABASE(id).DATAscanrate;
    %VMAXLIST=DATABASE(id).(regionvmaxlist);
    %VMAX=[VMAX;mean(VMAXLIST(:,3))];
    
    %determine maximum upstroke velocity for each trace
    UPSTROKERANGE=[floor(0.2*tracedatascanrate):ceil(1.0*tracedatascanrate)];
    VMAXLIST=[];
    for i=1:size(TRACEDATA,2)
        UPSTROKEVELOCITY=diff(TRACEDATA(UPSTROKERANGE,i))*tracedatascanrate;
        %get rid of spikes by median filtering
        MUPSTROKEVELOCITY=medfilt1(UPSTROKEVELOCITY-mean(UPSTROKEVELOCITY),5)+mean(UPSTROKEVELOCITY);    
        %average the 5 largest amplitudes
        SORTUPSTROKEAMPLITUDES=sort(MUPSTROKEVELOCITY,'descend');
        VMAXLIST=[VMAXLIST;mean(SORTUPSTROKEAMPLITUDES(1:5))];
    end
    
    VMAX=[VMAX;mean(VMAXLIST)];
    VMAXSTD=[VMAXSTD;std(VMAXLIST)];
    VMAXSEM=[VMAXSEM;std(VMAXLIST)/sqrt(size(VMAXLIST,1))];
    
    if isfield(DATABASE,(regionvmaxlistsm))==1 && isempty(DATABASE(id).(regionvmaxlistsm))==0
        VMAXSM=[VMAXSM;mean(VMAXLISTSM(:,3))];
        VMAXSTDSM=[VMAXSTDSM;std(VMAXLISTSM(:,3))];
        VMAXSEMSM=[VMAXSEMSM;std(VMAXLISTSM(:,3))/sqrt(size(VMAXLISTSM,1))];
    end
    
    
    end
    end
    waitbar(id/length(DATABASE));
end
%angle statistics
ANGSTAT=anglestat(ANGLISTN,22.5);
%ANG(i,:)=[BINCENTERS(i),meanvel,freq];
close (hdl)
fprintf('FINISHED!\n')
%% save results
newpath=uigetdir(loadpath,'Pick folder to save data');

if length(newpath)>2

    dataname=datafile(1:end-4);
    basefile=[num2str(dataname),'-',num2str(region)];

    if isempty(RATE)==0
        %save conduction data
        RATEOUTPUT=[ENTRY,RATE];
        ratepath=[newpath,'\',basefile,'-RATEVEL.txt'];
        save(ratepath,'RATEOUTPUT','-ascii');

        VELOUTPUT=[VEL,VELSTD];
        velpath=[newpath,'\',basefile,'-VEL.txt'];
        save(velpath,'VELOUTPUT','-ascii');
        
        DIVOUTPUT=[DIVPOSITIVE,DIVNEGATIVE,(DIVPOSITIVE+DIVNEGATIVE)];
        %DIV=[mean(DIVLIST(:,3))];
        %DIVPOSITIVE=[mean(DIVLISTPOSITIVE(:,3))]
        %DIVNEGATIVE=[mean(DIVLISTNEGATIVE(:,3))]
        divpath=[newpath,'\',basefile,'-DIV.txt'];
        save(divpath,'DIVOUTPUT','-ascii');
        
        ANGOUTPUT=[ANG,ANGSTD];
        angpath=[newpath,'\',basefile,'-MEANANG.txt'];
        anglistpath=[newpath,'\',basefile,'-ANGLISTN.txt'];
        angstatpath=[newpath,'\',basefile,'-ANGSTAT.txt'];
        save(angpath,'ANGOUTPUT','-ascii');
        save(anglistpath,'ANGLISTN','-ascii');
        save(angstatpath,'ANGSTAT','-ascii');
        
        POLOUTPUT=[POLROI,REPOLROI];%polarization times of rectangular ROI
        polpath=[newpath,'\',basefile,'-POLREPOLms.txt'];
        save(polpath,'POLOUTPUT','-ascii');
    end
    
    if isempty(RATEAPD)==0
        %save apd data
        RATEOUTPUT=[ENTRYAPD,RATEAPD];
        rateapdpath=[newpath,'\',basefile,'-RATEAPD.txt'];
        save(rateapdpath,'RATEOUTPUT','-ascii');

        APDOUTPUT=[APD,APDSTD,APDSEM];
        apdpath=[newpath,'\',basefile,'-APD.txt'];
        save(apdpath,'APDOUTPUT','-ascii');
        
        VMAXOUTPUT=[VMAX,VMAXSTD,VMAXSEM];
        vmaxpath=[newpath,'\',basefile,'-VMAX.txt'];
        save(vmaxpath,'VMAXOUTPUT','-ascii');
        
        if isfield(DATABASE,(regionvmaxlistsm))==1
            VMAXSMOUTPUT=[VMAXSM,VMAXSTDSM,VMAXSEMSM];
            vmaxsmpath=[newpath,'\',basefile,'-VMAXSM.txt'];
            save(vmaxsmpath,'VMAXSMOUTPUT','-ascii');
        end
    end
end
