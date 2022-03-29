%interpolate image stack

function [IDATA,iscanrate]=interpstack(DATA,scanrate,newscanrate)
%% run function 

%extract single pixel time series
if scanrate<newscanrate
    hdl = waitbar(0,['Interpolating Data']);
    %time axes
    TIME=[0:size(DATA,3)-1]'/scanrate;
    dt=1/newscanrate;
    NTIME=[TIME(1):dt:TIME(end)]';
    iscanrate=newscanrate;
    %initialize interpolated data matrix
    IDATA=zeros(size(DATA,1),size(DATA,2),length(NTIME));
else
    hdl = waitbar(0,['Downsampling Data']);
    n=floor(scanrate/newscanrate);
    iscanrate=scanrate/n;
    IDATA=zeros(size(DATA,1),size(DATA,2),ceil(size(DATA,3)/n));
end
for i=1:size(DATA,1)
    for j=1:size(DATA,2)
        PIXEL=squeeze(DATA(i,j,:));
        if scanrate<newscanrate
            IPIXEL=interp1(TIME,PIXEL,NTIME,'cubic');
            IDATA(i,j,:)=IPIXEL;
        else
            IPIXEL=downsample(PIXEL,n);
            IDATA(i,j,:)=IPIXEL;
        end
    end
    waitbar(i/size(DATA,1));
end
close(hdl)