function INVERTDATA=invertstack(DATA)
%% open Matlab Pool
%parallel processing tools
%%if matlabpool('size')==0
  %  h=msgbox('Activating MATLAB Pool','Parallel Processing Toolbox');
   % matlabpool open
    %close(h)
%end
%% read pixel data and invert
hdl = waitbar(0,['INVERTING IMAGE STACK']);
INVERTDATA=zeros(size(DATA));
for i=1:size(DATA,1)
    parfor j=1:size(DATA,2)
        PIXEL=squeeze(DATA(i,j,:));
        INVERTDATA(i,j,:)=(mean(PIXEL)-PIXEL)+mean(PIXEL);
    end
    waitbar(i/size(DATA,1));
end
close(hdl);