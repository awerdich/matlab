function hr=episodeheartrate(TRANSIENT,scanrate)
%% ALGORITHM TO CALCULATE HEART RATE FROM AP TRANSIENT
T=[1:size(TRANSIENT,1)]'/scanrate;
F=[0.01:0.001:5];
%define spectrum object
Hs=spectrum.periodogram;
%options for spectrum object
Hopts=msspectrumopts(Hs,(TRANSIENT-mean(TRANSIENT)));
%set custom options
Hopts.NormalizedFrequency=false;
Hopts.Fs=scanrate;
Hopts.SpectrumType='twosided';
Hopts.FreqPoints='User Defined';
Hopts.FrequencyVector=F;
%calculate mean-squared spectrum
%peaks reflect power in the signal at a given frequency band
Hmss=msspectrum(Hs,(TRANSIENT-mean(TRANSIENT)),Hopts);
%find maximum and normalize spectrum
SPEC=Hmss.Data;
[val,idx]=max(SPEC);NSPEC=SPEC/val;
 maxf=F(idx)*60;%heart rate in bpm