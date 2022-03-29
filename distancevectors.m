%estimate mean isochrone distance
function D=distancevectors
button=0;
D=[];
while button~=27
    if button~=27
        [AX,AY,p]=ginput(2);
        button=p(1);
        if button ~=27
            D=([D,[AX(2)-AX(1);AY(2)-AY(1)]])
        end
    end
end