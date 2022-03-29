%transfer angleunitvectors from one database to another database
%e.g. from CHAMBERS to OCIC to ensure that the same unitvectors are used in
%the angle calculations
%% load data files
if exist('loadpath','var')==0 || length(loadpath)<2
    [datafiles, loadpath] = uigetfile('*.mat', 'Pick a SOURCE DATABASE file');
    [datafilet, loadpath] = uigetfile('*.mat', 'Pick a TARGET DATABASE file');
end
load([loadpath,datafiles]);SDATABASE=DATABASE;
load([loadpath,datafilet]);TDATABASE=DATABASE;
clear DATABASE
%% find target filenames and replace angleunitvectors
PROCID=[];
for i=1:length(SDATABASE)
    sfile=SDATABASE(i).name;
    %find sfile in TDATABASE
    for j=1:length(TDATABASE)
        tfile=TDATABASE(j).name;
        if strcmp(sfile,tfile)==1
            %file found
            %replace unitvector in target database
            if isfield(TDATABASE,'ANGLEUNITVECTORXY') && sum(TDATABASE(j).ANGLEUNITVECTORXY==SDATABASE(i).ANGLEUNITVECTORXY)==2
                fprintf(['vector in file ',sfile,' in source database entry ',num2str(i),' already corrected \n']);
            else
                TDATABASE(j).ANGLEUNITVECTORXY=SDATABASE(i).ANGLEUNITVECTORXY;
                fprintf(['Replaced vector in file ',sfile,' in source database entry ',num2str(i),' .\n']);
                fprintf(['target file:',tfile,' source file:',sfile,' \n']);
                PROCID=[PROCID;i];
            end
        end
    end
end
%% remove old vector field in target database
if isfield(TDATABASE,'angleunitvectorxy')==1
    NTDATABASE=rmfield(TDATABASE,'angleunitvectorxy');
    TDATABASE=NTDATABASE;
end
%% save target database
file=input(['WRITE TO FILENAME:',num2str(datafilet),'>'],'s');
DATABASE=TDATABASE;
save(file,'DATABASE'); 