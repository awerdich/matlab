TRACE=load('CaG_testTL_Ic0_appended-right.txt','-ascii');
scanrate=125;
%% define fitting intervals
hdimage=figure('Name','FITTING INTERVALS');
plot(TRACE,'b');
XPOINTS=[];
%define baseline
button=1;
while button~=27
    [xi,yi,button]=ginput(1);
    if button~=27
        XPOINTS=[XPOINTS;round(xi)];
    end
end

%% mark fitting intervals
XDATA=[];INTERVALX=[];
figure(hdimage);hold on
for i=1:2:length(XPOINTS)-1
    INTERVALX=(XPOINTS(i):1:XPOINTS(i+1))';  
    XDATA=[XDATA;INTERVALX];
    plot(INTERVALX,TRACE(INTERVALX,:),'m')
end
%% data points for fitting
P=polyfit(XDATA,TRACE(XDATA),1);
FITXDATA=[1:length(TRACE)]';
FITYDATA=polyval(P,FITXDATA);
plot(FITXDATA,FITYDATA),'r';
FTRACE=TRACE-FITYDATA;
%% cut data for spectrum calculation
figure
plot(FTRACE);
[ax,ay]=ginput(2);
D=FTRACE(floor(min(ax)):floor(max(ax)));
%% spectrum
D=FTRACE(500:900);scanrate=125;
%set freqency vector [Hz]
F=[0:0.1:62.5]';
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
SPEC=Hmss.Data;      