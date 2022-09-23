function [CV,SELECTij]=cleanvfield(V,Vi,Vj,SELECTij,scale)
%V=[i,j,iabs,jabs,vel,ang]
figure
set(gca,'YDir','reverse');
hold on
for k=1:size(V,1)
    i=V(k,1);j=V(k,2);Vx=Vj(i,j);Vy=Vi(i,j);x=V(k,4);y=V(k,3);
    %plot selected velocity vectors
    quiver(x,y,Vx,Vy,scale,'LineWidth',1.5,'Color','k','MaxHeadSize',1);
    plot(x,y,'rs','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','g')
end
% mark already selected vectors
if exist('SELECTij')==1
    for s=1:size(SELECTij,1)
        m=SELECTij(s,1);n=SELECTij(s,2);
        quiver(n,m,Vj(m,n),Vi(m,n),scale,'LineWidth',1.5,'Color','r','MaxHeadSize',1);
    end
end
fprintf('select vectors to exclude in plot \n');

if exist('SELECTij')==0
    SELECTij=[];
end
button=1;
while button~=27
    [bx,by,button]=ginput(1);
    if isempty(button)==1
        button=27;
    end
    if button~=27
        %find coordinate
        i=round(by);j=round(bx);
        %check if selected vector is already on list
        onlist=0;%change only if selected coordinates found in list
        for k=1:size(SELECTij,1)
            if SELECTij(k,1)==i && SELECTij(k,2)==j
                %vector already selected
                onlist=k;
            end
        end
        if onlist==0
            SELECTij=[SELECTij;[i,j]];
            quiver(j,i,Vj(i,j),Vi(i,j),scale,'LineWidth',1.5,'Color','r','MaxHeadSize',1);
        elseif onlist>0
            m=SELECTij(onlist,1);n=SELECTij(onlist,2);
            quiver(n,m,Vj(m,n),Vi(m,n),scale,'LineWidth',1.5,'Color','k','MaxHeadSize',1);
            SELECTij(onlist,:)=[];
        end
    end
end
hold off
%find selected coordinates in VPLOT list and mark row
DELV=[];
for k=1:size(SELECTij,1)
    for l=1:size(V,1)
        if SELECTij(k,1)==V(l,1) && SELECTij(k,2)==V(l,2)
            DELV=[DELV;l];
        end
    end
end

%create new vectorlist an plot
CV=V;
CV(DELV,:)=[];
figure
set(gca,'YDir','reverse');
hold on
for k=1:size(CV,1)
    i=CV(k,1);j=CV(k,2);
    %plot selected velocity vectors
    quiver(j,i,Vj(i,j),Vi(i,j),scale,'LineWidth',1.5,'Color','k','MaxHeadSize',1);
    plot(j,i,'rs','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','g')
end
hold off
end