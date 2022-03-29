%estimate propagation velocity manually

button=0;
V=[];D=[];vi=[];
while button~=27
    if button~=27
        [AX,AY,p]=ginput(2);
        button=p(1);
        if button ~=27
            D=[AX,AY];
            pixeldistance=norm(D(2,:)-D(1,:),2);
            distance=pixeldistance*160/72*1.0e-3;%pixeldistance (mm)
            vi=distance/timeinterval*1.0e3;%velocity [mm/s]
            V=[V;vi]%
        end
    end
end
%statistics
MEANV=mean(V)
STDV=std(V)/sqrt(length(V))