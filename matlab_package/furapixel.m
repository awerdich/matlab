%open non-ratiometric ascii data and separate wavelengths
function [W1,W2]=furapixel(fratio,openmessage)

[file path]=uigetfile('*.txt',openmessage);
A=load([path,file]);
%% separate wavelengths 
ratios=floor(length(A)/fratio);
transcans=fratio/2-1;%number of frames for wavelength transition

for i=1:ratios
    f1=fratio*i-transcans-1;%frame for first wavelength
    f2=fratio*i;%frame for second wavelength
    W1(i)=A(f1);
    W2(i)=A(f2);
end


