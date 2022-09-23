%this function processed the ANGLISTN [ANGLE,CV] list
function ANG=anglestat(ANGLISTN,binwidth)
%% angle binning
%ANGLIST(i,:)=[vectorangle(REFERENCEVECTORij,VELOCITYVECTORij),VECTORLIST(i,5)];
ANGLEINTERVAL=[0,360];
firstbincenter=0;
bincenter=firstbincenter;
BINCENTERS=[];
INT=[];
while bincenter<ANGLEINTERVAL(2)-binwidth/2
    if bincenter-binwidth/2<0
        %first interval
        INT=[INT;[360-binwidth/2,bincenter+binwidth/2]];
    elseif 360<bincenter+binwidth/2
        %last interval
        INT=[INT;[bincenter-binwidth/2,360]];
    else
        INT=[INT;[bincenter-binwidth/2,bincenter+binwidth/2]];
    end
    BINCENTERS=[BINCENTERS;bincenter];
    bincenter=bincenter+binwidth;
end

%% collect velocity data for each bin
ANG=zeros(length(BINCENTERS),3);
for i=1:size(INT,1)
    %find angles within the INT angle interval
    ANGLIST_INT=[];
    for j=1:size(ANGLISTN,1)
        %first interval
        if INT(i,2)<INT(i,1)
            if (INT(i,1)<ANGLISTN(j,1) && ANGLISTN(j,1)<=ANGLEINTERVAL(2)) || (ANGLEINTERVAL(1)<ANGLISTN(j,1) && ANGLISTN(j,1)<=INT(i,2))
                ANGLIST_INT=[ANGLIST_INT;ANGLISTN(j,:)];
            end
        else
            if INT(i,1)<ANGLISTN(j,1) && ANGLISTN(j,1)<=INT(i,2)
                ANGLIST_INT=[ANGLIST_INT;ANGLISTN(j,:)];
            end
        end
    end
    %statistics for angles found in the INT angle interval
    
    if isempty(ANGLIST_INT)==1
        %there could be no vectors for a given angle interval
        meanvel=0;
        freq=0;
    else
        meanvel=mean(ANGLIST_INT(:,2));
        freq=100*size(ANGLIST_INT,1)/size(ANGLISTN,1);
    end
    %add statistics to output matrix
    ANG(i,:)=[BINCENTERS(i),meanvel,freq];
end
%add 360 degrees bin center
BINCENTERS=[BINCENTERS;ANGLEINTERVAL(2)];
INT=[INT;INT(1,:)];
ANG=[ANG;[ANGLEINTERVAL(2),ANG(1,2:end)]];
end