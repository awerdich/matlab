%Fura-2 Ca calibration
%% recording parameters and important other information
newcalibration=0;
if exist('CALIBRATIONDATA')==0
    newcalibration=1;
    fprintf('New calibration.\n');
    caldate=input('Calibration date YYMMDD>');
    camrate=input('Camera frame rate (frames/s)>');
    fratio=input('Frames per ratio>');
    lam1=input('First wavelength (nm)>');
    ban1=input('Bandwith of first wavelength (nm)>');
    lam2=input('Second wavelength (nm)>');
    ban2=input('Bandwidth of second wavelength (nm)>');
    repeat=1;
else
    askrepeat=input('Add new data to database (1=YES)?:');
    if askrepeat==1
        repeat=1;
    else
        repeat=0;
    end
end 
%% read Tiff stacks
while repeat==1
    %% process new image stack for each concentration
    if exist('CALIBRATIONDATA')==0
        id=1;%number of calibrations
        lastfile=[];
    else
        id=length(CALIBRATIONDATA)+1;
        lastfile=CALIBRATIONDATA(id-1).file;
    end
    %open first image
    if exist('path')==1
        defaultpath=path;
    else
        defaultpath=['.'];
    end
    [firstfile path]=uigetfile('*.tif',['IMAGE or STACK. Last:',num2str(lastfile)],defaultpath);
    if isempty(firstfile)==0 && length(firstfile)>1
    pathfirstfile=[path firstfile];

    %get file info from the first file
    info=imfinfo(num2str([path,firstfile]));
    
    %decide if selected file is a single image or a stack of many images
    if length(info)==1
        %new image stack
        fprintf('Building new image stack.\n');
    
    %estimate last frame number and number of digits used to count frames
    estlastframenumber=length(dir([path '*.tif']))-1;%disregard dark frame
    digits=size(num2str(estlastframenumber),2);%estimate number of digits
    lastframenumber=estlastframenumber;
    
    %reset first frame filename
    fprintf('Re-setting first file number to 1\n')
    %get correct first number
    firstnumber=sprintf(['%0',num2str(digits),'.0f'],1);
    firstfile=[firstfile(1:end-4-digits),firstnumber,firstfile(end-3:end)];
    firstframenumber=str2num(firstfile(end-3-digits:end-4));
    
    %correct last frame number to obtain ratios for all frames in file
    if mod(lastframenumber,fratio)~=0
        lastframenumber=floor(lastframenumber/fratio);
    end  
    %% create filename and get folder for image stacks
    %decide the number of digits required
    %create filename for each image and open them
    stackfile=[firstfile(1:end-4-digits),'-',num2str(firstframenumber),'-',num2str(lastframenumber),'.tif'];
    if exist('stackpath')==0
        %ask folder
        [f,stackpath]=uiputfile('*.tif','SELECT FOLDER FOR IMAGE STACK',defaultpath);clear f;
    end
    
    %delete old file, otherwise new frames will be added
    if exist([stackpath,stackfile])>0;
        fprintf('delete old image stack\n');
        return
    end
    %% assemble new image stack
    clear RAWDATA NRAWDATA;
    RAWDATA=double(zeros(info(1).Width,info(1).Height,lastframenumber-firstframenumber+1)); %pre-allocate 3D array 
    hdl = waitbar(0,['ASSEMBLING NEW IMAGE STACK ',num2str(stackfile)]);
    for i=firstframenumber:lastframenumber
        %add leading zeros to i and convert framenumber into string
        filei=[firstfile(1:end-4-digits),sprintf(['%0',num2str(digits),'d'],i),'.tif'];
        pathfilei=[path filei];
        j=i-firstframenumber+1;%frame number in data file
        %read frame
        newframe=uint16(imread(pathfilei,'tif'));
        %add to RAWDATA
        framenumber=i-firstframenumber+1;
        [RAWDATA(:,:,framenumber)]=double(newframe);%read data as double precision numbers
        %write frame into stack
        imwrite(newframe,[stackpath,stackfile],'tif','writemode','append','compression','none');
        waitbar((i-firstframenumber+1)/(lastframenumber-firstframenumber+1));
    end;
    close(hdl);
    %load dark frame
    %assemble file name
    DARKFRAME=uint16(imread([path,firstfile(1:end-4-digits),'_dark.tif']));
    %save copy of dark frame in stackpath
    imwrite(DARKFRAME,[stackpath,firstfile(1:end-4-digits),'_dark.tif']);
    DARKFRAME=double(DARKFRAME);%convert dark frame to double precision
    elseif length(info)>1 %selected image is a multiframe image stack
    %% read existing image stack
    %[firstfile path]=uigetfile('*.tif','FIRST IMAGE or STACK');
    stackfile=firstfile;
    stackpath=path;
    clear RAWDATA NRAWDATA;
    RAWDATA=double(zeros(info(1).Width,info(1).Height,length(info))); %pre-allocate 3D array 
    %read images into DATA array
    hdl = waitbar(0,['READING STACK ',num2str(stackfile)]);
    for frame=1:length(info)
        [RAWDATA(:,:,frame)]=double(imread([stackpath,stackfile],frame));
        waitbar(frame/length(info));
    end
    close(hdl)
    %load dark frame
    %assemble file name
    digits=size(num2str(length(info)),2);%digits of last frame number
    darkfile=[stackfile(1:end-digits-7),'_dark.tif'];
    %[darkfile darkpath]=uigetfile('*.tif',[num2str(stackfile),' SELECT DARK FRAME'],stackpath);
    DARKFRAME=double(imread([stackpath,darkfile],'tif'));
    end
    %New stack created or existing stack loaded
    %% Normalize RAWDATA array by DARKFRAME
    %correct for gain differences between recordings
    %pixel gains can vary between recordings; 
    %DARKFRAME is the average of 10 frames before shutter opening
    %DARKFRAME IS AN OFFSET AND NOT A GAIN!!!!!
    NRAWDATA=zeros(size(RAWDATA));
    for i=1:size(RAWDATA,3)
        NRAWDATA(:,:,i)=RAWDATA(:,:,i)-DARKFRAME;
    end
    %% separate wavelengths
    transcans=fratio/2-1;%number of frames for wavelength transition
    numframes=size(NRAWDATA,3)/fratio;
    W1=zeros(size(NRAWDATA,1),size(NRAWDATA,2),numframes);%first wavelength
    W2=zeros(size(NRAWDATA,1),size(NRAWDATA,2),numframes);%second wavelength
    R=zeros(size(NRAWDATA,1),size(NRAWDATA,2),numframes);%ratio
    for i=1:numframes
        f1=fratio*i-transcans-1;%frame for first wavelength
        f2=fratio*i;%frame for second wavelength
        W1(:,:,i)=NRAWDATA(:,:,f1);
        W2(:,:,i)=NRAWDATA(:,:,f2);
        R(:,:,i)=W1(:,:,i)./W2(:,:,i);
    end
    %mean and standard deviation of fluorescence values for each pixel
    MEANW1=zeros(size(W1,1),size(W1,2));STDW1=zeros(size(W1,1),size(W1,2));
    MEANW2=zeros(size(W2,1),size(W2,2));STDW2=zeros(size(W2,1),size(W2,2));
    MEANR=zeros(size(R,1),size(R,2));STDR=zeros(size(R,1),size(R,2));
    hdl = waitbar(0,'CALCULATING MEAN FLUORESCENCE VALUES');
    for i=1:size(R,1)
        for j=1:size(R,2)
            PIXELW1=squeeze(W1(i,j,:));
            PIXELW2=squeeze(W2(i,j,:));
            PIXELR=squeeze(R(i,j,:));
            MEANW1(i,j)=mean(PIXELW1);STDW1(i,j)=std(PIXELW1);
            MEANW2(i,j)=mean(PIXELW2);STDW2(i,j)=std(PIXELW2);
            MEANR(i,j)=mean(PIXELR);STDR(i,j)=std(PIXELR);
        end
        waitbar(i/size(R,1))
    end     
    close(hdl);
    %% save data in database
    %caldate=input('Calibration date YYMMDD>');
    %camrate=input('Camera frame rate (frames/s)>');
    %fratio=input('Frames per ratio>');
    %lam1=input('First wavelength (nm)>');
    %ban1=input('Bandwith of first wavelength (nm)>');
    %lam2=input('Second wavelength (nm)>');
    %ban2=input('Bandwidth of second wavelength (nm)>');
    stackfile
    concentration=input('Free Ca concentration (nmol/l)>');
    special=input('Measurement class <RETURN>=CONCENTRATION 1=BACKGROUND 2=0CA 3=SATURATED>');
    if isempty(special)==1
        special=0;
    end
    CALIBRATIONDATA(id).date=caldate;
    CALIBRATIONDATA(id).file=stackfile;
    CALIBRATIONDATA(id).path=stackpath;
    CALIBRATIONDATA(id).camrate=camrate;
    CALIBRATIONDATA(id).fratio=fratio;
    CALIBRATIONDATA(id).scanrate=camrate/fratio;
    CALIBRATIONDATA(id).lam1=lam1;
    CALIBRATIONDATA(id).ban1=ban1;
    CALIBRATIONDATA(id).lam2=lam2;
    CALIBRATIONDATA(id).ban2=ban2;
    CALIBRATIONDATA(id).concnanomol=concentration;
    CALIBRATIONDATA(id).special=special;
    CALIBRATIONDATA(id).specialhelp='0=CONCENTRATION 1=BACKGROUND 2=0CA 3=SATURATED';
    CALIBRATIONDATA(id).DATALAM1=W1;
    CALIBRATIONDATA(id).DATALAM2=W2;
    CALIBRATIONDATA(id).MEANLAM1=MEANW1;
    CALIBRATIONDATA(id).STDLAM1=STDW1;
    CALIBRATIONDATA(id).MEANLAM2=MEANW2;
    CALIBRATIONDATA(id).STDLAM2=STDW2;
    CALIBRATIONDATA(id).MEANR=MEANR;
    CALIBRATIONDATA(id).STDR=STDR;
    CALIBRATIONDATA(id).DARKFRAME=DARKFRAME;
    else
        repeat=0;%file selected empty. Exit main loop.
    end
end
%% Analyze CALIBRATIONDATA
%get concentrations and corrected ratios
%0=CONCENTRATION 1=BACKGROUND 2=0CA 3=SATURATED
%extract indices
%get wavelengths
lam1=CALIBRATIONDATA(1).lam1;
lam2=CALIBRATIONDATA(2).lam2;

CAIDX=[];BIDX=[];FREEIDX=[];SATIDX=[];
for i=1:length(CALIBRATIONDATA)
    if CALIBRATIONDATA(i).special==0
        CAIDX=[CAIDX;i];%Calcium data
    elseif CALIBRATIONDATA(i).special==1
        BIDX=[BIDX;i];%Background intensities
    elseif CALIBRATIONDATA(i).special==2
        FREEIDX=[FREEIDX;i];%Minimum intensity (free dye)
    elseif CALIBRATIONDATA(i).special==3
        SATIDX=[SATIDX;i];%Maximum intensity (saturated dye)
    end
end

%mean background intensities
if isempty(BIDX)==0
    BACKGROUND1=CALIBRATIONDATA(BIDX(1)).MEANLAM1;
    BACKGROUND2=CALIBRATIONDATA(BIDX(1)).MEANLAM2;
    if length(BIDX)>1
        for i=2:length(BIDX)
            BACKGROUND1=BACKGROUND1+CALIBRATIONDATA(BIDX(i)).MEANLAM1;
            BACKGROUND2=BACKGROUND2+CALIBRATIONDATA(BIDX(i)).MEANLAM2;
        end
        BACKGROUND1=BACKGROUND1./length(BIDX);
        BACKGROUND2=BACKGROUND2./length(BIDX);
    end
end

%minimum ratio Rmin
FREE1=CALIBRATIONDATA(FREEIDX(1)).MEANLAM1;
FREE2=CALIBRATIONDATA(FREEIDX(1)).MEANLAM2;
if length(FREEIDX)>1
    for i=2:length(FREEIDX)
        FREE1=FREE1+CALIBRATIONDATA(FREEIDX(i)).MEANLAM1;
        FREE2=FREE2+CALIBRATIONDATA(FREEIDX(i)).MEANLAM2;
    end
    FREE1=FREE1./length(FREEIDX);
    FREE2=FREE2./length(FREEIDX);
end
%CALCLUMATE MIN RATIO
if lam1<lam2
    RMIN=(FREE1-BACKGROUND1)./(FREE2-BACKGROUND2);
else
    RMIN=(FREE2-BACKGROUND2)./(FREE1-BACKGROUND1);
end

%maximum ratio Rmax
SAT1=CALIBRATIONDATA(SATIDX(1)).MEANLAM1;
SAT2=CALIBRATIONDATA(SATIDX(1)).MEANLAM2;
if length(SATIDX)>1
    for i=2:length(SATIDX)
        SAT1=SAT1+CALIBRATIONDATA(SATIDX(i)).MEANLAM1;
        SAT2=SAT2+CALIBRATIONDATA(SATIDX(i)).MEANLAM2;
    end
    SAT1=SAT1./length(SATIDX);
    SAT2=SAT2./length(SATIDX);
end
%CALCIULATE MAX RATIO
if lam1<lam2
    RMAX=(SAT1-BACKGROUND1)./(SAT2-BACKGROUND2);
else
    RMAX=(SAT2-BACKGROUND2)./(SAT1-BACKGROUND1);
end

%F0=F2MAX/F2MIN
if lam1<lam2
    F0=(FREE2-BACKGROUND2)./(SAT2-BACKGROUND2);
else
    F0=(FREE1-BACKGROUND1)./(SAT1-BACKGROUND1);
end

%Ca data
%sort CAIDX according concentrations
C=CALIBRATIONDATA(CAIDX(1)).concnanomol;
for i=2:length(CAIDX)
    C=[C;CALIBRATIONDATA(CAIDX(i)).concnanomol];
end
[SORTC,IDX]=sort(C);
CAIDX=CAIDX(IDX);%idices in CALIBRATIONDATA, sorted by Ca Concentrations

%calculate mean ratios
i=1;
CA=[];%final concentrations
R=[];%corresponding ratios

while i<=length(CAIDX)
    %get concentration
    IDX=CAIDX(i);%current index in CALIBRATIONDATA
    c=CALIBRATIONDATA(IDX).concnanomol;%corresponding concentration
    CA=[CA;c];%concentration
    %check for other entries with the same concentration
    %previous sorting ensures that entries with the same concentration
    %follow each other
    j=1;
    while (i+j)<=length(CAIDX) && CALIBRATIONDATA(CAIDX(i+j)).concnanomol==c
        IDX=[IDX;CAIDX(i+j)];
        j=j+1;
    end
    %calculate mean ratio for concentration c
    F1=CALIBRATIONDATA(IDX(1)).MEANLAM1;
    F2=CALIBRATIONDATA(IDX(1)).MEANLAM2;
    if lam1<lam2
        R0=(F1-BACKGROUND1)./(F2-BACKGROUND2);
    else
        R0=(F2-BACKGROUND2)./(F1-BACKGROUND1);
    end
    if length(IDX)>1
        for k=2:length(IDX)
            F1=CALIBRATIONDATA(IDX(k)).MEANLAM1;
            F2=CALIBRATIONDATA(IDX(k)).MEANLAM2;
            if lam1<lam2
                R0=R0+((F1-BACKGROUND1)./(F2-BACKGROUND2));
            else
                R0=R0+((F2-BACKGROUND2)./(F1-BACKGROUND1));
            end
        end
       R0=R0/length(IDX);
    end
    R(:,:,length(CA))=R0;
    i=i+j;
end
%% fit concentrations and ratios to estimate Kd

XDATA=log10(CA*1.0e-9);
fitpoints=1e4;
fraction=(XDATA(end)-XDATA(1))/fitpoints;
XFIT=[XDATA(1):fraction:XDATA(end)-fraction];

hdl = waitbar(0,'FITTING');
for i=1:size(R,1)
    for j=1:size(R,2)
        for k=1:length(CA)
            YDATA(k)=log10(F0(i,j)*(R(i,j,k)-RMIN(i,j))/(RMAX(i,j)-R(i,j,k)));
            %MODEL y=x+b <=> b=y-x;
            B(k)=YDATA(k)-XDATA(k);
        end
        b=mean(B);KD(i,j)=10^(-b)*1.0e9;
        %determine x-intercept graphically
        YFIT=XFIT+b;
        k=1;while YFIT(k)<0;k=k+1;end;xintercept=(XFIT(k-1)+XFIT(k))/2;KD(i,j)=(10^xintercept)*1.0e9;
    end
    waitbar(i/size(R,1));
end
close(hdl)
%% save fit results
calpath=uigetdir('Folder to save data files.');calpath=[calpath,'\'];
%save calibration data
filebase=['F2CAL',num2str(CALIBRATIONDATA(1).date)];
if newcalibration==1
    filedata=[filebase,'DATA'];
    save([calpath,filedata],'CALIBRATIONDATA');
end
fileresults=[filebase,'RESULTS'];
save([calpath,fileresults],'F0','KD','RMIN','RMAX','BACKGROUND1','BACKGROUND2');
    
    
