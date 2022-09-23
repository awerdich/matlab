%build tiff stack from single files
%set delfile=1 to delete origianl tif files
%get filename and path

%open first image
if exist('imagepath')==0 || length(imagepath)<2
    imagepath=[];firstfile=[];
end

[firstfile imagepath]=uigetfile({'*.tsm';'*.da'},['LAST FILE:',num2str(firstfile)],imagepath);
pathfirstfile=[imagepath firstfile];

%% file type
[filepath, name, ext] = fileparts(firstfile);
ext = ext(2:end);

if strcmp(ext,'da') == 1
    DATAST=buildstackfunction_da(imagepath,firstfile);
    loadstack=0;
elseif strcmp(ext,'tsm') == 1
    DATAST=buildstackfunction_tsm(imagepath,firstfile);
    loadstack=0;
else
    loadstack=1;
end
if loadstack==1
%open tif stack
%estimate last frame number and number of digits used to count frames
estlastframenumber=length(dir([imagepath '*.tif']))-1;%disregard dark frame
digits=size(num2str(estlastframenumber),2);%estimate number of digits

%ask if this is ratiometric data
if exist('askratio')==1 && askratio==1
    ratio=1;
else
    askratio=input('Is this ratiometric data? (1=YES)>');
    if askratio==1
        ratio=1;
    else
        ratio=0;
    end
end
    
if ratio==1
    camrate=500;fprintf(['Camera frame rate (frames/s)=',num2str(camrate),'\n'])  
    fratio=4;fprintf(['Frames per ratio=',num2str(fratio),'\n']);
    lam1=380;fprintf(['First wavelength (nm)=',num2str(lam1),'\n']);
    ban1=20;fprintf(['Bandwith of first wavelength (nm)=',num2str(ban1),'\n']);
    lam2=340;fprintf(['Second wavelength (nm)=',num2str(lam2),'\n']);
    ban2=20;fprintf(['Bandwidth of second wavelength (nm)=',num2str(ban2),'\n']);
    
    askparms=input('Change fura-2 recording parameters (1=yes, RETURN=no)>');
    if isempty(askparms)==0
        camrate=input('Camera frame rate (frames/s)>');
        fratio=input('Frames per ratio>');
        lam1=input('First wavelength (nm)>');
        ban1=input('Bandwith of first wavelength (nm)>');
        lam2=input('Second wavelength (nm)>');
        ban2=input('Bandwidth of second wavelength (nm)>');
    end
    fprintf('Re-setting first file number to 1\n')
    
    
    %get correct first number
    firstnumber=sprintf(['%0',num2str(digits),'.0f'],1);
    firstfile=[firstfile(1:end-4-digits),firstnumber,firstfile(end-3:end)];
end

%get info from the first file
info=imfinfo(num2str([imagepath,firstfile]));
firstframenumber=str2num(firstfile(end-3-digits:end-4));

%ask for last frame number suggesting estimate
inputlastframenumber=input(['lastframenumber [',num2str(estlastframenumber),']>']);
if isempty(inputlastframenumber)==0,lastframenumber=inputlastframenumber;else lastframenumber=estlastframenumber;end;

%correct last frame number to use ratios for all frames in file
if ratio==1 && mod(lastframenumber,fratio)~=0
    lastframenumber=floor(lastframenumber/fratio);
end
    

%read images into data matrix
n=lastframenumber-firstframenumber+1;
%decide the number of digits required
%create filename for each image and open them
%create filename for tif stack
if exist('stackpath')==0 || length(stackpath)<2
    stackpath=uigetdir(imagepath,'FOLDER FOR IMAGE STACK'); 
    stackpath=[stackpath,'\'];
end
stackfile=[firstfile(1:end-4-digits),'-',num2str(firstframenumber),'-',num2str(lastframenumber),'.tif'];

%delete old file, otherwise new frames will be added
if exist([stackpath,stackfile])>0;
        delete([stackpath,stackfile]);
        fprintf('delete file\n');
        return
end
pathstackfile=[stackpath,stackfile];
%% assemble stack
RAWDATA=double(zeros(info(1).Width,info(1).Height,lastframenumber-firstframenumber+1)); %pre-allocate 3D array 
hdl = waitbar(0,['LOADING IMAGE STACK ',num2str(stackfile)]);
for i=firstframenumber:lastframenumber
    %add leading zeros to i and convert framenumber into string
    filei=[firstfile(1:end-4-digits),sprintf(['%0',num2str(digits),'d'],i),'.tif'];
    pathfilei=[imagepath filei];
    j=i-firstframenumber+1;%frame number in data file
    %read frame
    newframe=uint16(imread(pathfilei,'tif'));
    %add to RAWDATA
    framenumber=i-firstframenumber+1;
    [RAWDATA(:,:,framenumber)]=double(newframe);
    %write frame into stack
    imwrite(newframe,pathstackfile,'tif','writemode','append','compression','none');
waitbar((i-firstframenumber+1)/(lastframenumber-firstframenumber+1));
end;
close(hdl);
%% load and save dark frame
DARKFRAME=uint16(imread([imagepath,firstfile(1:end-4-digits),'_dark.tif']));
%save copy of dark frame in stackpath
imwrite(DARKFRAME,[stackpath,firstfile(1:end-4-digits),'_dark.tif']);
DARKFRAME=double(DARKFRAME);%convert dark frame to double precision
%save fura parameters
if ratio==1
    file=[stackfile(1:end-4),'-F2PARMS.mat'];
    save([stackpath,file],'firstframenumber','lastframenumber','camrate','fratio','lam1','ban1','lam2','ban2');
end
end
