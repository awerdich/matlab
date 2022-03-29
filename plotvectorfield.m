%plots vectorfield from vectorcomponent matrices VX(i,j), VY(i,j)
function plotvectorfield(VX,VY)
%% plot vectorfield
vectorscale=5;
vectorwidth=1;
vectorcolor=[0 0 0];
%vector coordinates
X=[1:size(VX,2)];
Y=[1:size(VY,1)];
[XM,YM]=meshgrid(X,Y);
% plot vectorfield
figure
hold on
for i=1:2:length(Y)
    for j=1:2:length(X)
        quiver(X(j),Y(i),VX(i,j),VY(i,j),vectorscale,'LineWidth',vectorwidth,'Color',vectorcolor,'MaxHeadSize',vectorwidth);
    end
end
hold off
set(gca,'XLim',[1,length(X)],'YLim',[1,length(Y)],'YDir','reverse');