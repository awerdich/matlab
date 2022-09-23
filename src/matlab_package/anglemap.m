%function ANGMAP=anglemap(VECTORLIST,ANGLIST,VPLOTBIN)
id=1;
%this function needs the following data
%can be obtained by running getdatabase for a single sample
%VECTORLIST=[VECTORLIST;[i,j,vi,vj,norm([vi,vj],2)]];
%ANGLIST(i,:)=[vectorangle(REFERENCEVECTORij,VELOCITYVECTORij),VECTORLIST(i,5)];

%interpolated signal matrix
SIGNALFRAMEINTERP_0=DATABASE(id).v_SIGNALFRAMEINTERP0;
SIG=SIGNALFRAMEINTERP_0;
%angle map with all angles
ANGMAP=zeros(size(SIG));
%reduced angle map with 0<=angle<=180
ANGMAP180=zeros(size(SIG));
    for k=1:length(ANGLIST)
        i=VECTORLIST(k,1);j=VECTORLIST(k,2);
        angle=ANGLIST(k,1);
        if SIG(i,j)==1
            ANGMAP(i,j)=angle;
            if angle>180
                ANGMAP180(i,j)=360-angle;
            else
                ANGMAP180(i,j)=angle;
            end
        end
    end
%down sample data
[Xq,Yq]=meshgrid([1:2:size(SIG,2)],[1:2:size(SIG,1)]);
[X,Y]=meshgrid([1:size(SIG,2)],[1:size(SIG,1)]);
ANG=interp2(X,Y,ANGMAP,Xq,Yq);
ANG180=interp2(X,Y,ANGMAP180,Xq,Yq);
%prepare matrices for imscatter script
ASIGNALPIXELS=zeros(size(ANG));
ASIGNALPIXELS(ANG>0)=1;
[AS,map]=gray2ind(mat2gray(ASIGNALPIXELS));
RGBFBR=ind2rgb(AS,gray(2));
%add data for showmatrix
stackfile=DATABASE(1).name;
