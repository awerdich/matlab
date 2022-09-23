%% Equalize unit vectors in two databases

%opend databases
if exist('loadpath','var')==0 || length(loadpath)<2;loadpath=[];end

fprintf('Load database file with reference unit vectors (SOURCE DATABASE)\n');
[datafile1, loadpath] = uigetfile('*.mat', 'Pick a SOURCE DATABASE file',loadpath);
load([loadpath,datafile1]);SDATABASE=DATABASE;clear DATABASE;

fprintf('Load database file with vectors to overwrite (TARGET DATABASE)\n');
[datafile2, loadpath] = uigetfile('*.mat', 'Pick a TARGET DATABASE file',loadpath);
load([loadpath,datafile2]);TDATABASE=DATABASE;clear DATABASE;


%% examine target databases for matching entries in source database
REPLACEIDT=[];
TVECTOR=[];
for idt=1:length(TDATABASE)
    %look for matching filename in source database
    IDS=[];
    for ids=1:length(SDATABASE)
        %compare file names
        if strcmp(TDATABASE(idt).name,SDATABASE(ids).name)
            IDS=[IDS;ids];
        end
    end
    if isempty(IDS)==0 && length(IDS)==1
        %replace unit vector in target database
        TDATABASE(idt).ANGLEUNITVECTORXY=SDATABASE(IDS(1)).ANGLEUNITVECTORXY;
        fprintf(['Replaced vector in file:',SDATABASE(IDS(1)).name,' in database:',datafile2,'.\n']);
        REPLACEIDT=[REPLACEIDT;idt];
    elseif isempty(IDS)==0 && length(IDS)>1
        %warn that there are two database entries with the same filename
        fprintf(['File:',TDATABASE(idt).name,' has multiple entries in database:',datafile1,' IDs:',num2str(IDS),' .\n']);
        fprintf(['Skipped file:',TDATABASE(idt).name,'. \n']);
    end
    %collect final vector
    TVECTOR=[TVECTOR,TDATABASE(idt).ANGLEUNITVECTORXY];
end
%% summarize results
fprintf(['Replaced ',num2str(length(REPLACEIDT)),'/',num2str(length(TDATABASE)),' vectors.\n']);
datafile2a=[datafile2(1:end-4),'-EQ','.mat'];
DATABASE=TDATABASE;save(datafile2a,'DATABASE');
fprintf(['Saved target database ',datafile2a,' in current matlab folder.\n']);
TVECTOR
