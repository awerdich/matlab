%%Extract Values


if exist('N', 'var') == 1
   clear N
end

if exist('N', 'var') == 0
    
    %APD
    %varNames = {'NAME', 'PATH', 'V_APD', 'VOC_APD', 'VIC_APD', 'AV_APD', 'A_APD', 'AOC_APD', 'AIC_APD'};
    %varTypes = {'char', 'char', 'double','double','double','double','double','double','double'};
    %CV
    varNames = {'NAME', 'PATH', 'V_VEL', 'VOC_VEL', 'VIC_VEL', 'AV_VEL', 'A_VEL', 'AOC_VEL', 'AIC_VEL'};
    varTypes = {'char', 'char', 'double','double','double','double','double','double','double'};

    %CA
    %varNames = {'NAME', 'PATH', 'V VEL', 'VOC VEL', 'VIC VEL', 'AV VEL', 'A VEL', 'AOC VEL', 'AIC VEL'};

    
end
sz = [size(DATABASE,2) 9];
N = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
T = table('Size',[1 9],'VariableTypes',varTypes,'VariableNames',varNames);

        for i = 1:sz
            
            T.NAME = cellstr(DATABASE(i).name);
            T.PATH = cellstr(DATABASE(i).path);
            if isempty(DATABASE(i).v_vel) == 0
                T.V_VEL = DATABASE(i).v_vel;
            else
                T.V_VEL = 0;
            end
            
            if isempty(DATABASE(i).voc_vel) == 0
                T.VOC_VEL = DATABASE(i).voc_vel;
            else
                T.VOC_VEL = 0;
            end
            
            if isempty(DATABASE(i).vic_vel) == 0
                T.VIC_VEL = DATABASE(i).vic_vel;
            else
                T.VIC_VEL = 0;
            end
            
            if isempty(DATABASE(i).av_vel) == 0
               T.AV_VEL = DATABASE(i).av_vel;
            else
               T.AV_VEL = 0;
            end
            if isempty(DATABASE(i).a_vel) == 0
                T.A_VEL = DATABASE(i).a_vel;
            else
                T.A_VEL = 0;
            end
            
            if isempty(DATABASE(i).aoc_vel) == 0
                T.AOC_VEL = DATABASE(i).aoc_vel;
            else
                T.AOC_VEL = 0;
            end
            
            if isempty(DATABASE(i).aic_vel) == 0
                
                T.AIC_VEL = DATABASE(i).aic_vel;
            else
                T.AIC_VEL = 0;
            end
            
            N(i,:) = T;
        end
