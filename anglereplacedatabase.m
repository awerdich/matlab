%% replace unitvectors in older version of OC-IC database
% new database should only have a single unit vector per heart
% ANGLEUNITVECTORXY
VECTORS=[];
for i=1:length(DATABASE)
    if isempty(DATABASE(i).v_ANGLEUNITVECTORXY)==0
        VECTORS(i,:)=DATABASE(i).v_ANGLEUNITVECTORXY';
    else
        VECTORS(i,:)=DATABASE(i).a_ANGLEUNITVECTORXY';
    end
end
%% overwrite unitvectors
for id=1:length(DATABASE)
    fprintf(['FILENAME: ',num2str(DATABASE(id).name),'\n']);
    vline=input('corresponding line index in VECTORS array:')
    %calculate angle between old and new vectors
    fprintf('Old vector:\n'),DATABASE(id).angleunitvectorxy
    fprintf('New vector:\n'),VECTORS(vline,:)'
    OLDVECTOR=[DATABASE(id).angleunitvectorxy(2);DATABASE(id).angleunitvectorxy(1)];
    NEWVECTOR=[VECTORS(vline,2);VECTORS(vline,1)];
    angle360=vectorangle(OLDVECTOR,NEWVECTOR);
    fprintf(['angle between old and new vector:',num2str(angle360),' degrees \n']);
    input('Replace vector [RETURN] \n');
    DATABASE(id).ANGLEUNITVECTORXY=VECTORS(vline,:)';
end
%% double check filenames and vectors
for id=1:length(DATABASE)
    DATABASE(id).name
    DATABASE(id).ANGLEUNITVECTORXY
    VECTORS(id,:)'
    pause
end
%% remove angleunitvectorfields
%region=input('PIXELDATA region (sa, a, aic, aoc, av, v, vic, voc, o):','s');
if isfield(DATABASE,'angleunitvectorxy')==1
    DATABASE=rmfield(DATABASE,'angleunitvectorxy');
end
if isfield(DATABASE,'sa_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'sa_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'a_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'a_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'aic_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'aic_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'aoc_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'aoc_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'av_ANGLEUNTIVECTORXY')==1
    DATABASE=rmfield(DATABASE,'av_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'v_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'v_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'voc_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'voc_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'vic_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'vic_ANGLEUNITVECTORXY');
end
if isfield(DATABASE,'o_ANGLEUNITVECTORXY')==1
    DATABASE=rmfield(DATABASE,'o_ANGLEUNITVECTORXY');
end
