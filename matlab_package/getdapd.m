%determine APD gradient from DATABASE
%% load data file
if exist('loadpath','var')==0 || length(loadpath)<2
    [datafile, loadpath] = uigetfile('*.mat', 'Pick a DATABASE file');
end
load([loadpath,datafile]);
%% extract APDs from atrium and ventricle
APDGRADA=[];APDGRADV=[];APDGRADAV=[];%gradients
VMAXGRADA=[];VMAXGRADV=[];VMAXGRADAV=[];
RATE=[];HEARTID=[];
MEANAPDA=[];MEANAPDV=[];MEANAPDAV=[];
MEANVMAXA=[];MEANVMAXV=[];MEANVMAXAV=[];

for id=1:length(DATABASE)
%% process single data set
    %DATABASE filters
    if DATABASE(id).sinusloc==0 && DATABASE(id).age>48 
        
        %ATRIAL APDs
        APDRANGEA=[];APDA=[];
        VMAXRANGEA=[];VMAXA=[];
        
        if isempty(DATABASE(id).aic_APDLIST_ms)==0
            APD=DATABASE(id).aic_APDLIST_ms(:,3);
            VMAX=DATABASE(id).aic_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEA=[APDRANGEA;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),1]];
            VMAXRANGEA=[VMAXRANGEA;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),1]];
            APDA=[APDA;mean(APD)];
            VMAXA=[VMAXA;mean(VMAX)];
        end
    
        if isempty(DATABASE(id).a_APDLIST_ms)==0
            APD=DATABASE(id).a_APDLIST_ms(:,3);
            VMAX=DATABASE(id).a_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEA=[APDRANGEA;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),2]];
            VMAXRANGEA=[VMAXRANGEA;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),2]];
            APDA=[APDA;mean(APD)];
            VMAXA=[VMAXA;mean(VMAX)];
        end
    
        if isempty(DATABASE(id).aoc_APDLIST_ms)==0
            APD=DATABASE(id).aoc_APDLIST_ms(:,3);
            VMAX=DATABASE(id).aoc_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEA=[APDRANGEA;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),3]];
            VMAXRANGEA=[VMAXRANGEA;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),3]];
            APDA=[APDA;mean(APD)];
            VMAXA=[VMAXA;mean(VMAX)];
        end
   
        %VENTRICULAR APDs
        APDRANGEV=[];APDV=[];
        APDRANGEAV=[];APDAV=[];
        VMAXRANGEV=[];VMAXV=[];
        
        if isempty(DATABASE(id).voc_APDLIST_ms)==0
            APD=DATABASE(id).voc_APDLIST_ms(:,3);
            VMAX=DATABASE(id).voc_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEV=[APDRANGEV;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),3]];
            VMAXRANGEV=[VMAXRANGEV;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),3]];
            APDV=[APDV;mean(APD)];
            VMAXV=[VMAXV;mean(VMAX)];
        end
    
        if isempty(DATABASE(id).v_APDLIST_ms)==0
            APD=DATABASE(id).v_APDLIST_ms(:,3);
            VMAX=DATABASE(id).v_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEV=[APDRANGEV;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),2]];
            VMAXRANGEV=[VMAXRANGEV;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),2]];
            APDV=[APDV;mean(APD)];
            VMAXV=[VMAXV;mean(VMAX)];
        end
    
        if isempty(DATABASE(id).vic_APDLIST_ms)==0
            APD=DATABASE(id).vic_APDLIST_ms(:,3);
            VMAX=DATABASE(id).vic_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEV=[APDRANGEV;[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11)),1]];
            VMAXRANGEV=[VMAXRANGEV;[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11)),1]];
            APDV=[APDV;mean(APD)];
            VMAXV=[VMAXV;mean(VMAX)];
        end
    
        if isempty(DATABASE(id).av_APDLIST_ms)==0
            APD=DATABASE(id).av_APDLIST_ms(:,3);
            VMAX=DATABASE(id).vic_VMAXLISTSM_s(:,3);
            SORTAPD=sortrows(APD);
            SORTVMAX=sortrows(VMAX);
            %determine range
            APDRANGEAV=[mean(SORTAPD(11:20)),mean(SORTAPD(end-20:end-11))];
            VMAXRANGEAV=[mean(SORTVMAX(11:20)),mean(SORTVMAX(end-20:end-11))];
        end
        
        %Atrial data
        if isempty(APDRANGEA)==0
            APDGRADA=[APDGRADA;max(APDA)-min(APDA)];
            VMAXGRADA=[VMAXGRADA;max(VMAXA)-min(VMAXA)];
            MEANAPDA=[MEANAPDA;mean(APDA)];
            MEANVMAXA=[MEANVMAXA;mean(VMAXA)];
        end
        
        %Ventricular data
        if isempty(APDRANGEV)==0
            APDGRADV=[APDGRADV;max(APDV)-min(APDV)];
            VMAXGRADV=[VMAXGRADV;max(VMAXV)-min(VMAXV)];
            MEANAPDV=[MEANAPDV;mean(APDV)];
            MEANVMAXV=[MEANVMAXV;mean(VMAXV)];
        end
        
        %AV data
        APDGRADAV=[APDGRADAV;max(APDRANGEAV)-min(APDRANGEAV)];
        VMAXGRADAV=[VMAXGRADAV;max(VMAXRANGEAV)-min(VMAXRANGEAV)];
        MEANAPDAV=[MEANAPDAV;DATABASE(id).av_meanapd_ms20];
        MEANVMAXAV=[MEANVMAXAV;DATABASE(id).av_meanvmax_s];
        
        RATE=[RATE;DATABASE(id).rate];
        HEARTID=[HEARTID;id];
    end
end
save MEANAPDA.txt MEANAPDA -ascii -tabs
save MEANAPDV.txt MEANAPDV -ascii -tabs
save MEANVMAXA.txt MEANVMAXA -ascii -tabs
save MEANVMAXV.txt MEANVMAXV -ascii -tabs
save MEANAPDAV.txt MEANAPDAV -ascii -tabs
save MEANVMAXAV.txt MEANVMAXAV -ascii -tabs