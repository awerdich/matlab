%% GET F2 DATA
%% load database file
if exist('loadpath','var')==0 || length(loadpath)<2
    [datafile, loadpath] = uigetfile('*.mat', 'Pick a DATABASE file');
end
load([loadpath,datafile]);
%% extract data
fprintf(['First file in ',datafile,': ',FDATABASE(1).name,' .\n']);
region=input('ROI (sa, a, aic, aoc, av, v, vic, voc, o:','s');
clear regionv
if strcmp(region,'a')==1 || strcmp(region,'aic') || strcmp(region,'aoc') || strcmp(region,'sa')==1 ... 
|| strcmp(region,'av')==1 || strcmp(region,'v')==1 || strcmp(region,'o')==1 ...
|| strcmp(region,'vic')==1 || strcmp(region,'voc')==1 || strcmp(region,'vc')

    regionframedata=[region,'_DATA'];%ALL TRANSIENTS IN FRAME (i,j,:)
    regionmeantrace=[region,'_TRACE'];%MEAN FRAME TRANSIENT
    regionframex=[region,'_framex'];%frame width in pixels
    regionframey=[region,'_framey'];%frame height
    regionframeumwidth=[region,'_frameumwidth'];%frameumwidth
    regionframeumheight=[region,'_frameumheight'];%frameumheight
    regionframecenter=[region,'_framecenter'];%frame center locatin
    regionarea=[region,'_area'];%estimate of active tissue area in frame (SIGNALS>0)
    regionlow=[region,'_diast'];%AVERAGE (MINIMUM OF EACH TRANSIENT)
    regionhigh=[region,'_syst'];%AVERAGE (MAXIMUM OF EACH TRANSIENT)
    regionamp=[region,'_amp'];%AVERAGE (AMPLITUDE OF EACH TRANSIENT)
    regionframeij=[region,'_frameij'];%pixel coordinates for each transient used to calculate average
    regiondur=[region,'_durms50'];%AVERAGE TRANSIENT DURATION (ms) 
    regionlowmatrix=[region,'_RLOWROI'];%diastolic ROI matrix
    regionhighmatrix=[region,'_RHIGHROI'];%systolic ROI matrix
    regiondurmatrix=[region,'_DUROI'];%duration ROI matrix
else
    fprintf('Wrong input!\n');
    return;
end

DATA=[];TRANS=[];DATE=[];
DATACHAMBER=[];
RLOWLISTN=[];RHIGHLISTN=[];RAMPLISTN=[];DURLISTN=[];
for i=1:length(FDATABASE)
        %DATA [DIA,REL,DUR50]
        DATA=[DATA;[FDATABASE(i).(regionlow),FDATABASE(i).(regionamp),FDATABASE(i).(regiondur)]];
        DATE=[DATE;FDATABASE(i).date];
        %TRANSIENTS
        TRANS=[TRANS,FDATABASE(i).(regionmeantrace)];
        %get mean data for chambers
        if strcmp(region,'v')==1 && isfield(FDATABASE,'vic_diast')==1
            diast=mean([FDATABASE(i).vic_diast;FDATABASE(i).v_diast;FDATABASE(i).voc_diast]);
            amp=mean([FDATABASE(i).vic_amp;FDATABASE(i).v_amp;FDATABASE(i).voc_amp]);
            dur=mean([FDATABASE(i).vic_durms50;FDATABASE(i).v_durms50;FDATABASE(i).voc_durms50]);
            DATACHAMBER=[DATACHAMBER;[diast,amp,dur]];
        elseif strcmp(region,'a')==1 && isfield(FDATABASE,'aic_diast')==1
            diast=mean([FDATABASE(i).aic_diast;FDATABASE(i).a_diast;FDATABASE(i).aoc_diast]);
            amp=mean([FDATABASE(i).aic_amp;FDATABASE(i).a_amp;FDATABASE(i).aoc_amp]);
            dur=mean([FDATABASE(i).aic_durms50;FDATABASE(i).a_durms50;FDATABASE(i).aoc_durms50]);
            DATACHAMBER=[DATACHAMBER;[diast,amp,dur]];
        end
end

%save data
if exist('datapath','var')==0
    datapath=loadpath;
end
datapath=uigetdir(datapath,'SELECT FOLDER FOR DATA FILES'); 
datapath=[datapath,'\'];
save([num2str(datapath),datafile(1:end-4),'_',region,'_DATA.txt'],'DATA','-ascii');
save([num2str(datapath),datafile(1:end-4),'_',region,'_TRANS.txt'],'TRANS','-ascii');
save([num2str(datapath),datafile(1:end-4),'_',region,'_DATE.txt'],'DATE','-ascii');
if strcmp(region,'v')==1
    save([num2str(datapath),datafile(1:end-4),'_',region,'_CHAMBER.txt'],'DATACHAMBER','-ascii');
end