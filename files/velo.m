%determine velocities and angles for velocity components Vi(i,j),Vj(i,j)
function [VEL,ANG]=velo(Vi,Vj,calvfactor_i,calvfactor_j,UNITVECTOR)
VEL=zeros(size(Vi));ANG=zeros(size(Vi));
for i=1:size(Vi,1)
    for j=1:size(Vi,2)
        if Vi(i,j)~=0 && Vj(i,j)~=0
            VELVECTOR=[Vi(i,j)*calvfactor_i;Vj(i,j)*calvfactor_j];%velocity vector [mm/s] in matrix space
            v=norm(VELVECTOR,2);
            a=acos(dot(VELVECTOR,UNITVECTOR)/v)*180/pi;%angle in degrees
            if VELVECTOR(1)<0
                a=360-a;
            end
            VEL(i,j)=v;
            ANG(i,j)=a;
        end
    end
end
end