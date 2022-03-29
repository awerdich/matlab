%% eliminate signals outside the tissue area
opengl software
%convert signal image to RGB
[SIG,map]=gray2ind(ASIGNALPIXELS,2);
RGBSIG=ind2rgb(SIG,gray(2));

ALPHASIG=ASIGNALPIXELS*0.25;
%change color of signals
for i=1:size(ASIGNALPIXELS,1)
    for j=1:size(ASIGNALPIXELS,2)
        if ASIGNALPIXELS(i,j)==1
            RGBSIG(i,j,:)=[1,0,0];
        end
    end
end


figure('name','SELECT PIXELS TO EXCLUDE FROM FURTHER ANALYSIS');

%estimate pixels
%RGBFBR = DATAST(1).RGBFBR;
image(RGBFBR);hold on
image(RGBSIG,'AlphaData',ALPHASIG)
%set(gca,'XLim',[1 size(NORMDATA,2)],'YLim',[1 size(NORMDATA,1)]);

if exist('SCATTER')~=1 || isempty(SCATTER)
    SCATTER=ones(size(ASIGNALPIXELS));
end

AUTOSCATTER=ones(size(SCATTER));
%mark SCATTER==0 in image
for i=1:size(SCATTER,1)
    for j=1:size(SCATTER,2)
        if SCATTER(i,j)==0
            rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
        end
    end
end 
%click pixels until ESC
%button=0;
%save copy of SCATTER ARRAY
if exist('SCATTER')==1
    BACKUPSCATTER=SCATTER;
else
    BACKUPSCATTER=[];
end

while button~=27 && button ~=127
    clear ax ay
    [ax,ay,button]=ginput(1);
    if button~=1
        ax=[];ay=[];%if no pixel was selected, clear coordinates
    end
    if button~=27
        if isempty(ax)==0 && (button~=28 || button~=29 || button~=30 || button ~=31)
            %discard decimals 
            i=round(ay);j=round(ax); 
            %add position to the list of pixels
            if SCATTER(i,j)==1
                %pixel not marked
                SCATTER(i,j)=0;
                rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
                AUTOSCATTER(i,j)=0;
            elseif SCATTER(i,j)==0
                %pixel already marked
                SCATTER(i,j)=1;
                rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','k','FaceColor',RGBFBR(i,j,:));
                AUTOSCATTER(i,j)=1;
            end
        end

        if button==28 || button==29 || button==30 || button==31
            SCATTERLIST=[];
            for i=1:size(AUTOSCATTER,1)
                for j=1:size(AUTOSCATTER,2)
                    if AUTOSCATTER(i,j)==0;
                        SCATTERLIST=[SCATTERLIST;[i,j]];
                    end
                end
            end
        end
        
        %FILL EDGES
        if button==28
            %left edge
            for k=1:size(SCATTERLIST,1)
                i=SCATTERLIST(k,1);j=SCATTERLIST(k,2);
                while ASIGNALPIXELS(i,j)==1 && j>1
                    AUTOSCATTER(i,j)=0;
                    j=j-1;
                end
            end
            %add AUTOSCATTER frames to image
            for i=1:size(AUTOSCATTER,1)
                for j=1:size(AUTOSCATTER,2)
                    if AUTOSCATTER(i,j)==0
                        rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
                    end
                end
            end 
            %transfer AUTOSCATTER to SCATTER and reset AUTOSCATTER
            SCATTER(AUTOSCATTER==0)=0;
            AUTOSCATTER=ones(size(SCATTER));
        end

        if button==29
            %right edge
            for k=1:size(SCATTERLIST,1)
                i=SCATTERLIST(k,1);j=SCATTERLIST(k,2);
                while ASIGNALPIXELS(i,j)==1 && j<size(AUTOSCATTER,2)
                    AUTOSCATTER(i,j)=0;
                    j=j+1;
                end
            end
            %add AUTOSCATTER frames to image
            for i=1:size(AUTOSCATTER,1)
                for j=1:size(AUTOSCATTER,2)
                    if AUTOSCATTER(i,j)==0
                        rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
                    end
                end
            end 
            %transfer AUTOSCATTER to SCATTER and reset AUTOSCATTER
            SCATTER(AUTOSCATTER==0)=0;
            AUTOSCATTER=ones(size(SCATTER));
        end
        
        if button==30
            %top edge
            for k=1:size(SCATTERLIST,1)
                i=SCATTERLIST(k,1);j=SCATTERLIST(k,2);
                while ASIGNALPIXELS(i,j)==1 && i>1
                    AUTOSCATTER(i,j)=0;
                    i=i-1;
                end
            end
            %add AUTOSCATTER frames to image
            for i=1:size(AUTOSCATTER,1)
                for j=1:size(AUTOSCATTER,2)
                    if AUTOSCATTER(i,j)==0
                        rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
                    end
                end
            end 
            %transfer AUTOSCATTER to SCATTER and reset AUTOSCATTER
            SCATTER(AUTOSCATTER==0)=0;
            AUTOSCATTER=ones(size(SCATTER));
        end

        if button==31
            %bottom
            for k=1:size(SCATTERLIST,1)
                i=SCATTERLIST(k,1);j=SCATTERLIST(k,2);
                while ASIGNALPIXELS(i,j)==1 && i<size(AUTOSCATTER,1)
                    AUTOSCATTER(i,j)=0;
                    i=i+1;
                end
            end
            %add AUTOSCATTER frames to image
            for i=1:size(AUTOSCATTER,1)
                for j=1:size(AUTOSCATTER,2)
                    if AUTOSCATTER(i,j)==0
                        rectangle('Position',[j-0.5,i-0.5,1,1],'LineWidth',0.5,'EdgeColor','r','FaceColor','w');
                    end
                end
            end 
            %transfer AUTOSCATTER to SCATTER and reset AUTOSCATTER
            SCATTER(AUTOSCATTER==0)=0;
            AUTOSCATTER=ones(size(SCATTER));
        end
    end
end
%% save SCATTER data
if button==27 
    if exist('stackpath')==1 && isempty(datapath)==0
        [fpath, fname, fext] = fileparts(stackfile);
        scatterfile=[stackfile(1:end-length(fext)),'SCATTER.mat'];
        if exist('datapath')==0
            [file datapath]=uiputfile('*.*','SELECT FOLDER');
        end
        pathscatterfile=[datapath,scatterfile];
        save(pathscatterfile,'SCATTER','FR','RGBFBR','datapath')
    else
        save('SCATTER.mat','SCATTER');
    end
elseif button==127 
    SCATTER=BACKUPSCATTER;
end
