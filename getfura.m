A_BASELINE=[];
A_AMP=[];
A_MAX=[];
A_DUR=[];
V_BASELINE=[];
V_AMP=[];
V_MAX=[];
V_DUR=[];

for id=1:length(FDATABASE)
    date=FDATABASE(id).date;
    
    if isempty(FDATABASE(id).a_diast)==1 
        %atrium is empty
        a_diast=0;a_amp=0;a_syst=0;a_durms50=0;
    else
        a_diast=FDATABASE(id).a_diast;
        a_amp=FDATABASE(id).a_amp;
        a_syst=FDATABASE(id).a_syst;
        a_durms50=FDATABASE(id).a_durms50;
    end
    A_BASELINE=[A_BASELINE;[date,a_diast]];
    A_AMP=[A_AMP;[date,a_amp]];
    A_MAX=[A_MAX;[date,a_syst]];
    A_DUR=[A_DUR;[date,a_durms50]];
    
    if isempty(FDATABASE(id).v_diast)==1
        %ventricle is empty
        v_diast=0;v_amp=0;v_syst=0;v_durms50=0;
    else
        v_diast=FDATABASE(id).v_diast;
        v_amp=FDATABASE(id).v_amp;
        v_syst=FDATABASE(id).v_syst;
        v_durms50=FDATABASE(id).v_durms50;
    end
    V_BASELINE=[V_BASELINE;[date,v_diast]];
    V_AMP=[V_AMP;[date,v_amp]];
    V_MAX=[V_MAX;[date,v_syst]];
    V_DUR=[V_DUR;[date,v_durms50]];
end

save('A_BASELINE.txt','A_BASELINE','-ascii');
save('A_AMP.txt','A_AMP','-ascii');
save('A_MAX.txt','A_MAX','-ascii');
save('A_DUR.txt','A_DUR','-ascii');
save('V_BASELINE.txt','V_BASELINE','-ascii');
save('V_AMP.txt','V_AMP','-ascii');
save('V_MAX.txt','V_MAX','-ascii');
save('V_DUR.txt','V_DUR','-ascii');
