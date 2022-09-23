%Fura-2 Ca calibration
%% recording parameters and important other information
if exist('CDATA')==0 || isempty(CDATA)==1
fprintf('New calibration.\n');
caldate=input('Calibration date YYMMDD>');
camrate=input('Camera frame rate (frames/s)>');
fratio=input('Frames per ratio>');
lam1=input('First wavelength (nm)>');
ban1=input('Bandwith of first wavelength (nm)>');
lam2=input('Second wavelength (nm)>');
ban2=input('Bandwidth of second wavelength (nm)>');
%% Enter file names and concentrations
CDATA=[];i=1;%calibration data base
if exist('laststackfile')==0;laststackfile=[];laststackpath=[];end
while exist('stackpath')==0 || length(stackpath)>1
    [stackfile,stackpath]=uigetfile('*.tif',['SELECT IMAGE STACK. LAST: ',num2str(laststackfile)],laststackpath);
    if length(stackpath)>1
        laststackfile=stackfile;laststackpath=stackpath;
        %concentration
        stackfile
        concentration=input('Concentration (nM)>');
        repeat=1;
        while repeat==1
            type=input('(0) Background (no dye) (1) Concentration (2) Saturation [1]>');
            if isempty(type)==1
                type=1;repeat=0;
            elseif 0<=type && type<=2
                repeat=0;
            end
        end
        %save file info
        CDATA(i).date=caldate;
        CDATA(i).camrate=camrate;
        CDATA(i).fratio=fratio;
        CDATA(i).lam1=lam1;
        CDATA(i).ban1=ban1;
        CDATA(i).lam2=lam2;
        CDATA(i).ban2=ban2;
        CDATA(i).FILE=stackfile;
        CDATA(i).PATH=stackpath;
        CDATA(i).conc=concentration;
        CDATA(i).type=type;
    end
    i=i+1;
end
end
%% separate wavlengths and calculate mean fluorescence
for i=1:length(CDATA)
    %sort wavelengths in stack
    [W1,W2]=furasortwavelengths(CDATA(i).PATH,CDATA(i).FILE,fratio);
    %At 500frames/s, wavelength order is reversed in stack
    %Restore correct wavlength assignment for 500 frames/s
    if camrate==500
        W3=W1;W1=W2;W2=W3;clear W3;
    end
    %calculate mean fluorescence values
    CDATA(i).MEANW1=zeros(size(W1,1),size(W1,2));
    CDATA(i).MEANW2=zeros(size(W2,1),size(W2,2));
    CDATA(i).W1=W1;
    CDATA(i).W2=W2;
    for j=1:size(W1,1)
        for k=1:size(W1,2)
            CDATA(i).MEANW1(j,k)=mean(squeeze(W1(j,k,:)));
            CDATA(i).MEANW2(j,k)=mean(squeeze(W2(j,k,:)));
        end
    end
end
%% BACKGROUND, RMIN, RMAX
%type=input('(0) Background (no dye) (1) Concentration (2) Saturation [1]>');
%BACKGROUND
clear B1 B2 BACKGROUND1 BACKGROUND2
j=0;
for i=1:length(CDATA)
    if CDATA(i).type==0
        j=j+1;
        B1(:,:,j)=CDATA(i).MEANW1;
        B2(:,:,j)=CDATA(i).MEANW2;
    end
end 
BACKGROUND1=mean(B1,3);
BACKGROUND2=mean(B2,3);

%MINIMUM FLUORESCENCE
clear MIN1 MIN2 RMIN BMIN1 BMIN2
j=0;
for i=1:length(CDATA)
    if CDATA(i).type==1 && CDATA(i).conc==0
        j=j+1;
        MIN1(:,:,j)=CDATA(i).MEANW1;
        MIN2(:,:,j)=CDATA(i).MEANW2;
        %background-corrected minimum fluorescence
        BMIN1(:,:,j)=MIN1(:,:,j)-BACKGROUND1;
        BMIN2(:,:,j)=MIN2(:,:,j)-BACKGROUND2;
    end
end

%MAXIMUM FLUORESCENCE
clear MAX1 MAX2 RMAX BMAX1 BMAX2
j=0;
IDX=[];
for i=1:length(CDATA)
    if CDATA(i).type==2
        IDX=[IDX;i];
        j=j+1;
        MAX1(:,:,j)=CDATA(i).MEANW1;
        MAX2(:,:,j)=CDATA(i).MEANW2;
        %background-corrected maximum fluorescence
        BMAX1(:,:,j)=MAX1(:,:,j)-BACKGROUND1;
        BMAX2(:,:,j)=MAX2(:,:,j)-BACKGROUND2;
    end
end




















%% find negative values
NEG=[];
for i=1:size(BMAX1,1)
    for j=1:size(BMAX1,2)
        for k=1:size(BMAX1,3)
            if BMAX1(i,j,k)<BMIN1(i,j,k)
                NEG=[NEG;[i,j,k,BMAX1(i,j,k)]];
            end
        end
    end
end
%% RMIN RMAX
%both wavelengths belong together in each measurement
%AVERAGE RATIOS, NOT FLUORESCENCE VALUES!!!
R1=zeros(size(MAX1));
R2=zeros(size(MAX1));

RMAX=zeros(size(MAX1,1),size(MAX1,2));
RMIN=zeros(size(RMAX));
R0=zeros(size(RMAX));

if lam1<lam2
    for j=1:size(R,3)
        %minimum ratios
        R1(:,:,j)=BMIN1(:,:,j)./BMIN2(:,:,j);
        %maximum ratios
        R2(:,:,j)=BMAX1(:,:,j)./BMAX2(:,:,j);
    end
    %alpha=Sf2/Sb2
    R0=mean(BMIN2,3)./mean(BMAX2,3);
else
    for j=1:size(R,3)
        R1(:,:,j)=BMIN2(:,:,j)./BMIN1(:,:,j);
        R2(:,:,j)=BMAX2(:,:,j)./BMAX1(:,:,j);
    end
    R0=mean(BMIN1,3)./mean(BMAX1,3);
end
RMIN=mean(R1,3);
RMAX=mean(R2,3);
%% concentrations
%sorted list of concentrations
CALIST=[];
for i=1:length(CDATA)
    if CDATA(i).type==1 && 0<CDATA(i).conc
        CALIST=[CALIST;[i,CDATA(i).conc]];
    end
end
SORTCALIST=sortrows(CALIST,2);%[INDEX CDATA,CONCENTRATION]

%new database sorted by concentrations
SORTCALISTR=[];
%ratio frames
i=1;%first concentration
    SORTCALISTR(i).conc=CDATA(SORTCALIST(i,1)).conc;
    j=i;%count number of lines with same concentration
    CINDEX=SORTCALIST(j,:);%store first index

    %find lines with the same concentration
    while j+1<=length(SORTCALIST) && SORTCALIST(j,2)==SORTCALIST(j+1,2)
        j=j+1;
        CINDEX=[CINDEX;SORTCALIST(j,:)];
    end
    
    %calculate ratios for concentration i
    R=zeros(size(R0,1),size(R0,2),length(CINDEX));
    for j=1:size(R,3)
        idx=CINDEX(j,1);
        if lam1<lam2
            R(:,:,j)=(CDATA(idx).MEANW1-BACKGROUND1)./(CDATA(idx).MEANW2-BACKGROUND2);
        else
            R(:,:,j)=(CDATA(idx).MEANW2-BACKGROUND2)./(CDATA(idx).MEANW1-BACKGROUND1);
        end
    end
    
    %save mean ratio in new database
    SORTCALISTR(i).MEANR=mean








