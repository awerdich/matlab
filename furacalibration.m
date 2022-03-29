%Fura-2 Ca calibration
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
        CDATA(i).FILE=stackfile;
        CDATA(i).PATH=stackpath;
        CDATA(i).conc=concentration;
        CDATA(i).type=type;
    end
    i=i+1;
end
%% separate wavlengths and calculate mean fluorescence
for i=1:length(CDATA)
    %load fura parameters for each CDATA entry
    stackpath=num2str(CDATA(i).PATH);
    stackfile=num2str(CDATA(i).FILE);
    file=[stackfile(1:end-4),'-F2PARMS.mat'];
    load([stackpath,file],'firstframenumber','lastframenumber','camrate','fratio','lam1','ban1','lam2','ban2');
    %sort wavelengths in stack
    [W1,W2]=furasortwavelengths(stackpath,stackfile,fratio);
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
%% BACKGROUND MIN MAX DATA
%type=input('(0) Background (no dye) (1) Concentration (2) Saturation [1]>');
%BACKGROUND
clear B1 B2
BFILE=[];
j=0;
for i=1:length(CDATA)
    if CDATA(i).type==0
        j=j+1;
        BFILE(j).FILE=CDATA(i).FILE;
        B1(:,:,j)=CDATA(i).MEANW1;
        B2(:,:,j)=CDATA(i).MEANW2;
    end
end 

%CA FREE SOLUTION
clear MIN1 MIN2
MINFILE=[];
j=0;
for i=1:length(CDATA)
    if CDATA(i).type==1 && CDATA(i).conc==0
        j=j+1;
        MINFILE(j).FILE=CDATA(i).FILE;
        MIN1(:,:,j)=CDATA(i).MEANW1;
        MIN2(:,:,j)=CDATA(i).MEANW2;
    end
end

%CA SATURATED SOLUTION
clear MAX1 MAX2
MAXFILE=[];
j=0;
for i=1:length(CDATA)
    if CDATA(i).type==2
        j=j+1;
        MAXFILE(j).FILE=CDATA(i).FILE;
        MAX1(:,:,j)=CDATA(i).MEANW1;
        MAX2(:,:,j)=CDATA(i).MEANW2;
    end
end
%% MINIMUM AND MAXIMUM RATIOS
RMIN=zeros(size(MIN1,1),size(MIN1,2));
RMAX=zeros(size(RMIN));
F0=zeros(size(RMIN));
NEGRMAX2=[];
for i=1:size(MIN1,1)
    for j=1:size(MIN1,2)
        if lam1<lam2
            %380nm fluorescence at maximum ca (MAX2) can be near background
            %avoid MAX-B2 to be below background
            RMIN1=mean(squeeze(MIN1(i,j,:)))-mean(squeeze(B1(i,j,:)));
            RMIN2=mean(squeeze(MIN2(i,j,:)))-mean(squeeze(B2(i,j,:)));
            
            RMAX1=mean(squeeze(MAX1(i,j,:)))-mean(squeeze(B1(i,j,:)));
            RMAX2=mean(squeeze(MAX2(i,j,:)))-mean(squeeze(B2(i,j,:)));

            F1=mean(squeeze(MIN2(i,j,:)))-mean(squeeze(B2(i,j,:)));
            F2=mean(squeeze(MAX2(i,j,:)))-mean(squeeze(B2(i,j,:)));
            
            if RMAX2<0
               RMAX2=mean(squeeze(MAX2(i,j,:)))-min(squeeze(B2(i,j,:)));
               F2=mean(squeeze(MAX2(i,j,:)))-min(squeeze(B2(i,j,:)));
               NEGRMAX2=[NEGRMAX2;[i,j,RMAX2]];
            end
            
        else
            %380nm fluorescence at maximum ca (MAX1) can be near backround
            %avoid MAX1-B1 to be below background
            RMIN1=mean(squeeze(MIN2(i,j,:)))-mean(squeeze(B2(i,j,:)));
            RMIN2=mean(squeeze(MIN1(i,j,:)))-mean(squeeze(B1(i,j,:)));
            
            RMAX1=mean(squeeze(MAX2(i,j,:)))-mean(squeeze(B2(i,j,:)));
            RMAX2=mean(squeeze(MAX1(i,j,:)))-mean(squeeze(B1(i,j,:)));

            F1=mean(squeeze(MIN1(i,j,:)))-mean(squeeze(B1(i,j,:)));
            F2=mean(squeeze(MAX1(i,j,:)))-mean(squeeze(B1(i,j,:)));
            
            if RMAX2<0
                RMAX2=mean(squeeze(MAX1(i,j,:)))-min(squeeze(B1(i,j,:)));
                F2=mean(squeeze(MAX1(i,j,:)))-min(squeeze(B1(i,j,:)));
                NEGRMAX2=[NEGRMAX2;[i,j,RMAX2]];
            end
        end
        RMIN(i,j)=RMIN1/RMIN2;
        RMAX(i,j)=RMAX1/RMAX2;
        F0(i,j)=F1/F2;
    end
end
%replace pixels with negative ratios by filtered ratios
%spatial filter
FRMIN=wiener2(RMIN,[2 2]);
FRMAX=wiener2(RMAX,[2 2]);
FF0=wiener2(F0,[2 2]);
for i=1:size(NEGRMAX2,1)
    j=NEGRMAX2(i,1);k=NEGRMAX2(i,2);
    RMAX(j,k)=FRMAX(j,k);
    F0(j,k)=FF0(j,k);
end
%% CA Files

%measure dissociation constant
%collect indeces and concentrations
CALIST=[];
for i=1:length(CDATA)
    if CDATA(i).type==1 && 0<CDATA(i).conc
        CALIST=[CALIST;[i,CDATA(i).conc]];
    end
end

%continue only if there are calibration data
if isempty(CALIST)==0
    %sort CALIST
    SCALIST=sortrows(CALIST,2);

    %collect mean ratios
    CAFRAMES=[];
    i=1;%line in SCALIST
    cnum=0;%concentration number
    while i<size(SCALIST,1)    
    j=i;
    %find last line in SCALIST with same concentration
    while j+1<=size(SCALIST,1) && SCALIST(i,2)==SCALIST(j+1,2)
       j=j+1;
    end
   
    %collect frames from i to j
    W1=zeros(size(F0,1),size(F0,2),j-i+1);
    W2=zeros(size(W1));
    cnum=cnum+1;
    for k=i:j
       fprintf(['[',num2str(cnum),'] ',...
           num2str(CDATA(SCALIST(k,1)).conc),' nM ',...
           'Processing data file:',num2str(CDATA(SCALIST(k,1)).FILE),' \n']);
       W1(:,:,k-i+1)=CDATA(SCALIST(k,1)).MEANW1;
       W2(:,:,k-i+1)=CDATA(SCALIST(k,1)).MEANW2;
   end
   
   %calculate ratios and save results in new database
   CAFRAMES(cnum).conc=SCALIST(i,2);
   for m=1:size(W1,1)
    for n=1:size(W1,2)
        if lam1<lam2
            R1=mean(squeeze(W1(m,n,:)))-mean(squeeze(B1(m,n,:)));
            R2=mean(squeeze(W2(m,n,:)))-mean(squeeze(B2(m,n,:)));
        else
            R1=mean(squeeze(W2(m,n,:)))-mean(squeeze(B2(m,n,:)));
            R2=mean(squeeze(W1(m,n,:)))-mean(squeeze(B1(m,n,:)));
        end
        CAFRAMES(cnum).R(m,n)=R1./R2;
    end
   end
   i=j+1;
    end
    
    %analyze results
    %concentrations
    CA=[];B=[];
    for i=1:length(CAFRAMES)
        CA(i)=CAFRAMES(i).conc;
    end
    XDATA=log10(CA*1.0e-9);

    for i=1:size(CAFRAMES(1).R,1)
        for j=1:size(CAFRAMES(1).R,2)
            for k=1:length(CA)
                YDATA(k)=log10((F0(i,j))*(CAFRAMES(k).R(i,j)-RMIN(i,j))/(RMAX(i,j)-CAFRAMES(k).R(i,j)));
                %MODEL y=x+b <=> b=y-x;
                B(k)=YDATA(k)-XDATA(k);
            end
            b=mean(B);KD(i,j)=10^(-b)*1.0e9;
        end
    end
else
    KD=ones(size(RMIN));
end

%% save results
%save calibration data

BACKGROUND1=mean(B1,3);
BACKGROUND2=mean(B2,3);
file='CALIBRATION.mat';
save([stackpath,file],'F0','KD','RMIN','RMAX','BACKGROUND1','BACKGROUND2');