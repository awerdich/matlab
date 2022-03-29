%% startparameter
M=POL;%activation matrix
minfitpixels=8;%minimum number of pixels used in fit
maxframe=100;
maxrmse=0.4;%return only successful fits with rmse<minrmse
pixelcalfactor=16/7.201613;%REDSHIRT 80x80 [um/pixel]
calvfactor=scanrate*1.0e-3*pixelcalfactor;%[frame*mm/(pixel*s)][pixel/frame*frame/s*um/pixel*1000mm/um=mm/s]
estvel=20/calvfactor;%estimated mean velocity
%% remove timing signals marked by SCATTER and set FITPIXELS
if exist('SCATTER')==1
    for i=1:size(SCATTER,1)
        M(SCATTER(i,1),SCATTER(i,2),:)=0;%remove pixel marked by imregist
    end
end
%% get pixel from activation matrix
hdl = waitbar(0,'Fitting polynomial surfaces');
FITPARAMETERLIST=[];FITLOCATIONLIST=[];
FITPARAMETERARRAY=zeros(size(M));
RMSE=zeros(size(M));
for i=1:size(M,1)
    for j=1:size(M,2)
        if M(i,j)>0
%% fit pixels on wavefront
            %FIT
            LOCATION=[i,j];
            [Vij,parameter,rmse,pdistance]=fitpixels(M,LOCATION,estvel,minfitpixels);
            if 0<rmse && rmse<maxrmse %fit was successful
                %save fit parameters
                FITPARAMETERLIST=[FITPARAMETERLIST;parameter];%save fitparameters
                FITLOCATIONLIST=[FITLOCATIONLIST;LOCATION];%save LOCATION
                FITPARAMETERARRAY(i,j)=length(FITPARAMETERLIST);%locate p in FITPARAMETERLIST
                RMSE(i,j)=rmse;
            end
%% repeat fit for next pixel
        end
    end
waitbar(i/size(M,1));
end
close(hdl);
%% create new array of fitted activation times
maxpixeldistance=4;
CM=M;
hdl = waitbar(0,'Trying to correct activation times');
for i=1:size(M,1)
    for j=1:size(M,2)
        if RMSE(i,j)>0 && RMSE(i,j)<maxrmse %if fit error was small, replace activation time
            p=FITPARAMETERLIST(FITPARAMETERARRAY(i,j));
            CM(i,j)=polyvaln(p,[j,i]);%image coordinates
        elseif (RMSE(i,j)==0 || RMSE(i,j)>maxrmse) && M(i,j)>0
            %if fit error was larger than maxrmse, search for close pixel with small error
            ERROR=[];PAR=[];
            for k=1:size(FITLOCATIONLIST,1)
                m=FITLOCATIONLIST(k,1);n=FITLOCATIONLIST(k,2);
                %find entries in LOWERROR that are close to (i,j)
                if norm([j-n;i-m],2)<maxpixeldistance
                    ERROR=[ERROR;[RMSE(m,n)]];
                    PAR=[PAR;FITPARAMETERLIST(FITPARAMETERARRAY(m,n))];
                end
            end
            %from the collected pixels choose the one w/ the smallest error
            if isempty(ERROR)==0
                [val,idx]=min(ERROR);
                p=PAR(idx);
                CM(i,j)=polyvaln(p,[j,i]);
            end
        end
    end
    waitbar(i/size(M,1))
end
close(hdl)
%% save fit data
[actfile actpath]=uiputfile('*.mat','Pick folder');
actfile=[stackfile(1:end-4),'-ATIME.mat'];

fprintf(['saving activation time data \n'])
[actpath,actfile]
save([actpath,actfile],'POL','REPOL','M','CM');