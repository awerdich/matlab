%Showme(database)
[baseName, folder] = uigetfile({'*.da';'*.tsm'});
fullFileName = fullfile(folder, baseName)
load(strcat(folder,baseName(1:end-3),'-N.mat'));
load(strcat(folder,baseName(1:end-3),'-APD80.mat'));
load(strcat(folder,baseName(1:end-5),'SCATTER.mat'));
askshow = input('Showmatrix(1) or Showcontours(2)\n');
bo = 0;

if askshow == 1;
    %while bo == 0;
        showmatrix;
        %ax = gca;
        %bo = ax.XLim;
    %end;
else askshow == 2;  
    %while ax.XLim == 0;
       showcontours;
       %ax = gca;
       %bo = ax.XLim;
    %end;
end

askdescript = input('DESCRIPTION\n', 's');
saveas(gcf, strcat(folder,baseName,askdescript,'.png'));
clear all;    
    
