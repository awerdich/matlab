%combine multiple data arrays with offset correction

file=1;
counter=1;
%load first episode
[file path]=uigetfile('*.mat',['OPEN -N FILE EPISODE #1',num2str(counter)]);
load([path,file]);


while file~=0 
    % save previous arrays
    NORMDATA1=NORMDATA;
    if isempty(STIM)==1
        STIM1=zeros(size(NORMDATA,3),1);
    else
        STIM1=STIM;
    end
    FITINTERVALS1=FITINTERVALS;
    ASIGNALPIXELS1=ASIGNALPIXELS;
    EPISODEFRAME1=EPISODEFRAME;
    % load file
    counter=counter+1;
    [file path]=uigetfile('*.mat',['OPEN -N FILE EPISODE #',num2str(counter)]);
    if file~=0
        load([path,file],'NORMDATA','ASIGNALPIXELS','STIM','FITINTERVALS','EPISODEFRAME')
        EPISODEFRAME=[EPISODEFRAME1;EPISODEFRAME];
        FITINTERVALS2=FITINTERVALS;
        ASIGNALPIXELS(ASIGNALPIXELS1==0)=0;
        NORMDATA2=zeros(size(NORMDATA1,1),size(NORMDATA1,2),size(NORMDATA1,3)+size(NORMDATA,3));
        STIM=[STIM1;STIM];
        %STIM=[STIM1;STIM-(mean(STIM(FITINTERVALS2(1):FITINTERVALS2(2)))-mean(STIM1(FITINTERVALS1(end-2):FITINTERVALS1(end))))];
        %combine NORMDATA
        hdl = waitbar(0,['APPENDING DATA']);
        for i=1:size(ASIGNALPIXELS,1)
            for j=1:size(ASIGNALPIXELS,2)
                if ASIGNALPIXELS(i,j)>0
                    
                    PIXEL1=squeeze(NORMDATA1(i,j,:));
                    PIXEL2=squeeze(NORMDATA(i,j,:));
                    %PIXEL=[PIXEL1;PIXEL2-(PIXEL2(1)-PIXEL1(end))];
                    PIXEL=[PIXEL1;PIXEL2-(mean(PIXEL2(FITINTERVALS2(1):FITINTERVALS(2)))-mean(PIXEL1(FITINTERVALS1(end-1):FITINTERVALS1(end))))];
                    NORMDATA2(i,j,:)=PIXEL;
                    
                
                end
            end
            waitbar(i/size(ASIGNALPIXELS,1))
        end
        close(hdl)
        NORMDATA=NORMDATA2;
    end
end
%% save new data file
cfile=[stackfile(1:end-4),'-COMBINED',num2str(min(EPISODEFRAME(:))),'-',num2str(max(EPISODEFRAME(:))),num2str('-N.mat')];
if exist('datapath')==0
    [file datapath]=uiputfile('*.*','SELECT FOLDER');
end
stackpath=datapath;

save([datapath,cfile],'STIM','NORMDATA','ASIGNALPIXELS','scanrate','stackpath','datapath','cfile','stackfile','firstframenumber','lastframenumber','FBR','RGBFBR');