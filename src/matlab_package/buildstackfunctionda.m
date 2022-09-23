function [DATAST]=buildstackfunctionda(path,file)
%% get path and filename
if isempty(path)==1 || isempty(file)==1
    [file path]=uigetfile('*.da','Open .da file');
end
%% open file
fid=fopen([path,file],'r');
p='int16';
%% read header information
frewind(fid);%re-set file pointer to beginning of file 
    fread(fid,1,p);%POSITION 1
trials=fread(fid,1,p);%POSITION 2
    fread(fid,2,p);%POSITION 4
frames=double(fread(fid,1,p));%POSITION 5
    fread(fid,379,p);%POSITION 384
columns=double(fread(fid,1,p));%POSITION 385
rows=double(fread(fid,1,p));%POSITION 386
    fread(fid,2,p);%POSITION 388
scanrate=1/double(fread(fid,1,p))*1.0e6;%POSITION 389
    fread(fid,2,p);%POSITION 391
bncratio=double(fread(fid,1,p));%POSITION 392
    double(fread(fid,1,p));%POSITION 393
lam1=double(fread(fid,1,p));%POSITION 394
lam2=double(fread(fid,1,p));%POSITION 395
inslit1=double(fread(fid,1,p));%POSITION 396
inslit2=double(fread(fid,1,p));%POSITION 397
outslit1=double(fread(fid,1,p));%POSITION 398
outslit2=double(fread(fid,1,p));%POSITION 399
ratiorate=double(fread(fid,1,p));%POSITION 400
%% warning for multiple trials
if trials>1
    uiwait(msgbox(['Found ',num2str(trials),' trials in data file. OK to proceed.'],'Multiple trials'));
    asktrials=input(['Number of trials to export (<',num2str(trials),'). <ENTER> export all:']);
    if isempty(asktrials)==1
        extrials=trials;
    else
        if asktrials>trials
            extrials=trials;
            fprintf(['Exporting all',num2str(extrials),' trials.\n']);
        else
            extrials=asktrials;
            fprintf(['Exporting ',num2str(extrials),' trials.\n']);
        end
    end
else
    extrials=1;
end
%% load OPTICAL data
frewind(fid);%re-set file pointer to beginning of file
fread(fid,2560,p);%set file pointer to beginning of data
%POSITION 2560
for id=1:extrials
%hdl = waitbar(0,['Loading CCD data for trial:',num2str(id),'. Please wait.']);
%the bar takes just too long to process
if lam1>0
    fprintf(['Loading ratiometric CCD data for trial ',num2str(id),'. Please wait.\n']);
else
    fprintf(['Loading CCD data for trial ',num2str(id),'. Please wait.\n']);
end
RAWDATA=zeros(rows,columns,frames);
for k=1:(columns*rows)
    i=k/columns;%colum position of pixel k
    if i-floor(i)==0;%remainder of x/columns
        i=floor(i);
    else
        i=floor(i)+1;
    end
    j=k-columns*(i-1);%row position of pixel k
    RAWDATA(i,j,1:frames)=fread(fid,frames,p);%read pixel k
    %waitbar(k/(columns*rows));
end
%close(hdl)
%POSITION 2560+COLUMNS*ROWS*FRAMES
%% load BNC data
echan=8;%number of BNC channels
BNCDATA=zeros(frames*bncratio,echan);
for k=1:echan
    BNCDATA(:,k)=fread(fid,frames*bncratio,p);%read BNC data for channel k
end
%POSITION 2560 + COLUMNS*ROWS*FRAMES + ECHAN*FRAMES*BNCRATIO 
%% load DARKFRAME
DARKFRAME=zeros(rows,columns);
for k=1:(columns*rows)
    i=k/columns;%colum position of pixel k
    if i-floor(i)==0;%i at the end of a row
        i=floor(i);
    else
        i=floor(i)+1;
    end
    j=k-columns*(i-1);%row position of pixel k
    DARKFRAME(i,j)=fread(fid,1,p);%read pixel k
end
fread(fid,8,p);
%DARK FRAME IS STORED LAST (ROWS * COLUMNS + 8)
%% save data in structure DATAST
%project name
k=length(file);while strcmp(file(k),'.')==0 && k>1;k=k-1;end;filext=file(k+1:end);
DATAST(id).FILE=file;
DATAST(id).PATH=path;
DATAST(id).FILETYPE=filext;
DATAST(id).PROJECTNAME=file(1:k-1);
DATAST(id).trial=id;
DATAST(id).columns=columns;
DATAST(id).rows=rows;
DATAST(id).scanrate=scanrate;
DATAST(id).bncratio=bncratio;
DATAST(id).lam1=lam1;
DATAST(id).ban1=outslit1;
DATAST(id).lam2=lam2;
DATAST(id).ban2=outslit2;
DATAST(id).fratio=scanrate/ratiorate;
DATAST(id).RAWDATA=RAWDATA(1:80, 1:80, :);
DATAST(id).BNCDATA=BNCDATA;
DATAST(id).DARKFRAME=DARKFRAME(1:80, 1:80);
DATAST(1)
end
