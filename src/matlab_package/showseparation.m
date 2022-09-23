
%% show wavelength separation after findwav
line=90;%line in SCANINTERVALS to display
m=40;n=40;%matrix coordinates to dispolay

I=[min(SCANINTERVALS(line,:))-10:max(SCANINTERVALS(line,:))+10];
FLI=squeeze(RAWDATA(m,n,I));%corresponding fluroescence data

SI1=[SCANINTERVALS(line,1):SCANINTERVALS(line,2)];%first wavelength scans
FL1=squeeze(RAWDATA(m,n,SI1));%corresponding fluorescence data
SI2=[SCANINTERVALS(line,3):SCANINTERVALS(line,4)];%second wavelength scans
FL2=squeeze(RAWDATA(m,n,SI2));
%% plot WAVELENGTH DATA
wavfigure = figure('Name','WAVELENGTH SIGNAL');
plot(I,WAV(I),'b');hold on;plot(I,WAV(I),'m+');
plot(SI1,WAV(SI1),'ko');
plot(SI2,WAV(SI2),'ro');
%% plot READY DATA if available
if isempty(READY)==0
    readyfigure = figure('Name','READY SIGNAL');
    plot(I,READY(I),'b');hold on;plot(I,READY(I),'m+');
    plot(SI1,READY(SI1),'ko');
    plot(SI2,READY(SI2),'ro');
end
%% plot RAWDATA
fluorescencefigure = figure('Name','FLUORESCENCE SIGNAL');
plot(I,FLI,'b');hold on;plot(I,FLI,'m+');
plot(SI1,FL1,'ko');
plot(SI2,FL2,'ro');