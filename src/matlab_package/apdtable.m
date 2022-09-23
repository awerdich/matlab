%%Extract Values


if exist('N', 'var') == 1
   clear N
end

if exist('N', 'var') == 0
    
    %APD
    varNames = {'NAME', 'PATH', 'V_APD', 'VOC_APD', 'VIC_APD', 'AV_APD', 'A_APD', 'AOC_APD', 'AIC_APD'};
    varTypes = {'char', 'char', 'double','double','double','double','double','double','double'};
    %CV
    %v = 2; varNames = {'NAME', 'PATH', 'V VEL', 'VOC VEL', 'VIC VEL', 'AV VEL', 'A VEL', 'AOC VEL', 'AIC VEL'};
     %      varTypes = {'char', 'char', 'double','double','double','double','double','double','double'};

    %CA
    %varNames = {'NAME', 'PATH', 'V VEL', 'VOC VEL', 'VIC VEL', 'AV VEL', 'A VEL', 'AOC VEL', 'AIC VEL'};

    
end
sz = [size(DATABASE,2) 9];
N = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
T = table('Size',[1 9],'VariableTypes',varTypes,'VariableNames',varNames);

        for i = 1:sz
            
            T.NAME = cellstr(DATABASE(i).name);
            T.PATH = cellstr(DATABASE(i).path);
            if isempty(DATABASE(i).v_meanapd_ms20) == 0
                T.V_APD = DATABASE(i).v_meanapd_ms20;
            else
                T.V_APD = 0;
            end
            
%             if isempty(DATABASE(i).voc_meanapd_ms20) == 0
%                 T.VOC_APD = DATABASE(i).voc_meanapd_ms20;
%             else
%                 T.VOC_APD = 0;
%             end
            
%             if isempty(DATABASE(i).vic_meanapd_ms20) == 0
%                 T.VIC_APD = DATABASE(i).vic_meanapd_ms20;
%             else
%                 T.VIC_APD = 0;
%             end
%             
            if isempty(DATABASE(i).av_meanapd_ms20) == 0
                T.AV_APD = DATABASE(i).av_meanapd_ms20;
            else
                T.AV_APD = 0;
            end
            if isempty(DATABASE(i).a_meanapd_ms20) == 0
                T.A_APD = DATABASE(i).a_meanapd_ms20;
            else
                T.A_APD = 0;
            end
            
%             if isempty(DATABASE(i).aoc_meanapd_ms20) == 0
%                 T.AOC_APD = DATABASE(i).aoc_meanapd_ms20;
%             else
%                 T.AOC_APD = 0;
%             end
%             
%             if isempty(DATABASE(i).aic_meanapd_ms20) == 0
%                 
%                 T.AIC_APD = DATABASE(i).aic_meanapd_ms20;
%             else
%                 T.AIC_APD = 0;
%             end
            
            N(i,:) = T;
        end

    