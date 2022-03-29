function ANGLIST=getang(Vi,Vj,RMSE)
%% calculate velocity vector angles from old database file
ANGLIST=[];%list of velocity angles
UNITVECTOR=[0;1];
for i=1:size(Vi,1)
    for j=1:size(Vj,2)
        if 0<RMSE(i,j) && RMSE(i,j)<0.5  
        %absolute velocity [pixel/frame]    
        VELVECTOR=[Vi(i,j);Vj(i,j)];
        VELVECTOR=VELVECTOR/norm(VELVECTOR,2);%normalizes vector
        %angle
        a=acos(dot(VELVECTOR,UNITVECTOR))*180/pi;%angle in degrees
        if VELVECTOR(1)<0
            a=360-a;
        end
        ANGLIST=[ANGLIST;a];%matrix coordinates        
        end
    end
end
end