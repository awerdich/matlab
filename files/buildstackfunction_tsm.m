function [DATAST]=buildstackfunction_tsm(path,file)

% File parameters
userecords = 15; % Used records
maxrecord = 36; % Maximum number of records in header 
recordsize = 80; % length of each record (bytes)

% Get path and filename
if isempty(path)==1 || isempty(file)==1
    [file path]=uigetfile('*.tsm','Open .tsm file');
end

% Open file
fid=fopen([path,file],'r', 'l');

%% Get header
frewind(fid) % Set pointer to beginning of file
header = [];
for i=1:userecords
    s = fread(fid, [1, recordsize], '*char');
    header = [header; s];
end

%% Get values from the header

trials = 1;
lam1 = 380;
lam2 = 340;
inslit1 = 20;
inslit2 = 20;
outslit1 = 20;
outslit2 = 20;
ratiorate = 125;

[k, v] = getKeyVal(header(6,:));
frames = str2num(v);
[k, v] = getKeyVal(header(4, :));
rows = str2num(v);
[k, v] = getKeyVal(header(5, :));
columns = str2num(v);
[k, v] = getKeyVal(header(11, :));
scanrate = 1/str2num(v);

%% Load optical data
frewind(fid);
fseek(fid, 2880, 0); 

RAWDATA = zeros(rows,columns,frames);

for k = 1:frames
    
    f = fread(fid, [rows, columns], 'int16');
    RAWDATA(:,:,k) = f;
end

DARKFRAME = fread(fid, [rows, columns], 'int16');

%% Load BNC data

fid2=fopen([path,[file(1:end-4),'.tbn']],'r', 'l');
frewind(fid2)
echan = abs(fread(fid2, 1, 'int16'));
bncratio = fread(fid2, 1, 'int16');
scans = frames * bncratio; % Number of scans for each channel

BNCDATA = zeros(scans, echan);
for k = 1:echan
    BNCDATA(:, k) = fread(fid2, scans, 'double');
end

%% Put it all together
id = 1;% just one trial
DATAST(id).FILE=file;
DATAST(id).PATH=path;
DATAST(id).FILETYPE=file(end-2:end);
DATAST(id).PROJECTNAME=file(1:end-4);
DATAST(id).trial=trials;
DATAST(id).columns=columns;
DATAST(id).rows=rows;
DATAST(id).scanrate=scanrate;
DATAST(id).bncratio=bncratio;
DATAST(id).lam1=lam1;
DATAST(id).ban1=outslit1;
DATAST(id).lam2=lam2;
DATAST(id).ban2=outslit2;
DATAST(id).fratio=scanrate/ratiorate;
DATAST(id).RAWDATA=RAWDATA;
DATAST(id).BNCDATA=BNCDATA;
DATAST(id).DARKFRAME=DARKFRAME;
DATAST(1)

end

%% Local functions

function [key,val] = getKeyVal(rec)
    kv = strsplit(rec, '=');
    key = char(kv(1));
    key = key(find(~isspace(key)));
    if strcmp(key,'END') == 0
        val = char(kv(2));
        val = val(find(~isspace(val)));
    else
        val = [];
    end
end

