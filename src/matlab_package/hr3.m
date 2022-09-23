%% HERT RATES MEASUREMENT
%open data file
if exist('path')==0 || isempty('path')==1 || length(path)<2;
    path=[];
end
if exist('lastfile')==0
        lastfile=[];
end
selectfile=1;
DATA=[];i=1;
while selectfile==1
    if exist('lastpath')==1 && length(path)<2
        path=lastpath;
    end
    [file,path]=uigetfile('*.LOG',['open Metamorph .LOG files. Last file:',num2str(lastfile)],path);%get filename
    if length(file)>1
        lastpath=path;%save path
        lastfile=file;%save file
        %Read data
        D=importdata([path,file],',');
        %convert time column to numbers
        TSTRING=D.textdata(2:end,end-1);
        %get software version
        v=version('-release');
        vn=str2num(v(1:end-1));
        T=zeros(size(TSTRING));
        for j=1:length(T)
            %for new Matlab Version
            if vn>2011
                T(j)=str2num(TSTRING{j});
            else
                S=TSTRING{j};
                T(j)=str2num(S(2:end-1));
            end
        end
        DATA(i).TIME=T*1000;
        DATA(i).INTENSITY=D.data;
        i=i+1;
    else
        selectfile=0;
    end
end
%% Calculate dominant frequency for each dataset
hdl = waitbar(0,['CALCULATING RATES']);
RATE=[];
for j=1:length(DATA)
    TIME=DATA(j).TIME;
    IN=DATA(j).INTENSITY;
    %% establish equally spaced time axis by interpolation
    scanrate=100;%frames per second
    T=[TIME(1):1/scanrate*1.0e3:TIME(end)]';
    D=[interp1(TIME,IN,T,'cubic')];
    D=D-min(D);
    D=D/max(D);
    %% frequency analysis
    %set freqency vector [Hz]
    F=[0.01:0.001:5];
    %define spectrum object
    Hs=spectrum.periodogram;
    %options for spectrum object
    Hopts=msspectrumopts(Hs,(D-mean(D)));
    %set custom options
    Hopts.NormalizedFrequency=false;
    Hopts.Fs=scanrate;
    Hopts.SpectrumType='twosided';
    Hopts.FreqPoints='User Defined';
    Hopts.FrequencyVector=F;
    %calculate mean-squared spectrum
    %peaks reflect power in the signal at a given frequency band
    Hmss=msspectrum(Hs,(D-mean(D)),Hopts);
    %find maximum and normalize spectrum
    SPEC=Hmss.Data;
    [val,idx]=max(SPEC);NSPEC=SPEC/val;
    maxf=F(idx)*60;%heart rate in bpm
    %save frequency spectrum and rate in data file
    DATA(j).TIME=T;
    DATA(j).INTENSITY=D;
    DATA(j).SPECTRUM=NSPEC;
    DATA(j).F=F;
    DATA(j).ratebpm=maxf;
    RATE=[RATE;maxf];
    waitbar(j/length(DATA));
end
close(hdl);
%% display data
RATE
%% save data
sfile=[lastfile(1:end-4),'-DATA.mat'];
sfiletxt=[lastfile(1:end-4),'-RATESbpm.txt'];
save([lastpath,sfile],'DATA');
save([lastpath,sfiletxt],'RATE','-ascii')