%find signals in high frequency data with no baseline
%% initialize arrays and define parameters
DATA=FILTERDATA(:,:,1:250);
ASIGNALPIXELS=zeros(size(DATA,1),size(DATA,2));
minamp=110;
AMPLIST=[];
% calculate signal amplitudes
for i=1:size(DATA,1)
    for j=1:size(DATA,2)
% find signal amplitudes for single pixel
        PIXEL=squeeze(DATA(i,j,:));
        amp=max(PIXEL(150:end))-min(PIXEL(150:end));
        AMPLIST=[AMPLIST;amp];
        if amp>minamp
            ASIGNALPIXELS(i,j)=1;
        end
% start over with next pixel
    end
end
imagesc(ASIGNALPIXELS);colormap(gray)
%% normalize data
NDATA=double(zeros(size(DATA)));%normalized data array
%process ASIGNALPIXELS
for i=1:size(ASIGNALPIXELS,1)
    for j=1:size(ASIGNALPIXELS,2)
        if ASIGNALPIXELS(i,j)==1
            %extract pixel
            PIXEL=squeeze(DATA(i,j,:));
            %remove offset
            PIXEL=PIXEL-min(PIXEL);
            %normalize 
            PIXEL=PIXEL/max(PIXEL);
            %save in NORMDATA array
            NDATA(i,j,:)=PIXEL;
        end
    end
end
%% interpolate data
scanrate=125;%current scanrate
nscanrate=2000;%new scanrate
ndatalength=round((size(DATA,3)-1)*nscanrate/scanrate)+1;
NORMDATA=zeros(size(DATA,1),size(DATA,2),ndatalength);

hdl = waitbar(0,['interpolating DATA']);
for i=1:size(NDATA,1)
    for j=1:size(NDATA,2)
        %extract time series
        PIXEL=squeeze(NDATA(i,j,:));
        %create time axis
        xdata=(0:length(PIXEL)-1)/scanrate;
        nxdata=xdata(1):1/nscanrate:xdata(end);
        %interpolate PIXEL
        nydata=interp1(xdata,PIXEL,nxdata,'pchip');
        NORMDATA(i,j,:)=nydata';
    end
    waitbar(i/size(NDATA,1));
end
close(hdl)
            