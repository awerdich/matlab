%consistency check
VPLOTBIN=DATABASE(1).v_PLOTVELOCITYBIN;
CALVFRAMEj=DATABASE(1).v_VFRAMEj*DATABASE(1).scanrate*DATABASE(1).pixelcalfactorj*1.0e-3;
CALVFRAMEi=DATABASE(1).v_VFRAMEi*DATABASE(1).scanrate*DATABASE(1).pixelcalfactori*1.0e-3;
VLIST=[];
for i=1:size(VPLOTBIN,1)
        for j=1:size(VPLOTBIN,2)
            if VPLOTBIN(i,j)==1
                vx=CALVFRAMEj(i,j);
                vy=CALVFRAMEi(i,j);
                VLIST=[VLIST;norm([vx,vy],2)];
            end
        end
    end