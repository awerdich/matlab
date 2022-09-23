%save activation and repolarization time in database
%confirm that the correct heart is loaded
fprintf(['heart in database:',num2str(DATABASE(id).name),' \n']);
fprintf(['heart loaded:',num2str(stackfile),' \n']);
fprintf(['add field to database <RETURN>\n']);

%verify that DATABASE entry and loaded data belong to the same file
name=DATABASE(id).name(1:length(stackfile)-4);
ldata=stackfile(1:end-4);

if exist('stackfile','var')==1 && strcmp(name,ldata)==1
    correctfile=1;
    fprintf(['DATABASE file name:',num2str(DATABASE(id).name),'\n']);
    fprintf(['Data file name:',num2str(stackfile),'\n']);
    fprintf('DATABASE entry and loaded data confirmed!\n')
else
    correctfile=0;
    fprintf(['DATABASE file name:',num2str(DATABASE(id).name),'\n']);
    fprintf(['Data file name:',num2str(stackfile),'\n']);
    fprintf('Abort program. Check file names!\n');
end

if correctfile==1

    %SIGNAL MATRIX
    SIGNAL=ASIGNALPIXELS;
    SIGNAL(SCATTER==0)=0;
    SIGNAL(POL==0)=0;
    ROI=SIGNAL;
    
    %save POL and REPOL matrices
    DATABASE(id).POL=POL;
    DATABASE(id).REPOL=REPOL;
    DATABASE(id).SIGNAL=SIGNAL;
    
    %loop to get the ROI data
    
    for r=1:8
        %look up region
        if r==1;region='a';end
        if r==2;region='aic';end
        if r==3;region='aoc';end
        if r==4;region='av';end
        if r==5;region='v';end
        if r==6;region='vic';end
        if r==7;region='voc';end
        if r==8;region='o';end

        %define field names
        regionframex=[region,'_XFRAMEINTERP0'];%x-coordinates of frame region
        XFRAMEINTERP_0=DATABASE(id).(regionframex);
        regionframey=[region,'_YFRAMEINTERP0'];%y-coordinates of frame region
        YFRAMEINTERP_0=DATABASE(id).(regionframey);
        regionsignalframeinterp=[region,'_SIGNALFRAMEINTERP0'];%Binary image: Signals in frame
        SIGNALFRAMEINTERP_0=DATABASE(id).(regionsignalframeinterp);
        regionpol=[region,'_POL'];
        regionrepol=[region,'_REPOL'];
    
        %make sure ROI data is available
        if isfield(DATABASE,regionsignalframeinterp)==1 && isempty(DATABASE(id).(regionsignalframeinterp))==0
    
            %get the coordinates and initialize data ROIs
            X=[XFRAMEINTERP_0(1,1):XFRAMEINTERP_0(1,end)]';
            Y=[YFRAMEINTERP_0(1,1):YFRAMEINTERP_0(end,1)]';
            POL_ROI=zeros(length(Y),length(X));
            REPOL_ROI=zeros(size(POL_ROI));
    
            %get times
            for i=Y(1):Y(end)
                k=i-Y(1)+1;
                for j=X(1):X(end)
                    l=j-X(1)+1;
                    if SIGNAL(i,j)==1
                        ROI(i,j)=0;
                        POL_ROI(k,l)=POL(i,j);
                        REPOL_ROI(k,l)=REPOL(i,j);
                    end
                end
            end
    
            %save data in DATABASE
            DATABASE(id).(regionpol)=POL_ROI;
            DATABASE(id).(regionrepol)=REPOL_ROI;
        end
    end
end
DATABASE(id)
    
    
  