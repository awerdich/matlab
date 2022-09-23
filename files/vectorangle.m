function angle360=vectorangle(ANGLEUNITVECTORij,VELOCITYVECTORij)
%% determines the angle between angleuntivector and velocityvector [vi,vj]
%all input vectors in matrix space
%all calculations in image space
%convert angles to image space and normalize
N=[ANGLEUNITVECTORij(2);ANGLEUNITVECTORij(1)]/norm(ANGLEUNITVECTORij,2);
V=[VELOCITYVECTORij(2);VELOCITYVECTORij(1)]/norm(VELOCITYVECTORij,2);
E=[1;0];%unit vector in x-direction
beta=acos(dot(N,E));
%convention: 
%1. y-axis down
%2. angles are measured counterclockwise
if N(2)>0 
    beta=2*pi-beta;
end
alpha=acos(dot(V,E));
if V(2)>0
    alpha=2*pi-alpha;
end
gamma=alpha-beta;
%add 2pi if angle is negative
if gamma<0
    angle360=(gamma+2*pi)*180/pi;
else
    angle360=gamma*180/pi;
end
% %% plot velocity vectors
% VECTORLIST=[VECTORLIST;[i,j,vi,vj,norm([vi,vj],2)]];
% vectorwidth=2;
% allvectorfigure=figure;
% hold on
% %plot all velocity vectors in ROI
% for k=1:2:size(VECTORLIST,1)
%     x=VECTORLIST(k,2);y=VECTORLIST(k,1);vx=VECTORLIST(k,4);vy=VECTORLIST(k,3);
%     quiver(x,y,vx,vy,0.5,'LineWidth',vectorwidth,'Color','b','MaxHeadSize',vectorwidth);
% end
% %plot mean vector
% %MEANVECTORij=[sum(VECTORLIST(:,3)),sum(VECTORLIST(:,4))]/size(VECTORLIST,1);
% x=round(size(PLOTBIN,2)/2);y=round(size(PLOTBIN,1)/2);vx=MEANVECTORij(2);vy=MEANVECTORij(1);
% quiver(x,y,vx,vy,0.5,'LineWidth',vectorwidth,'Color','r','MaxHeadSize',vectorwidth);  
% %plot selected vector
% x=VECTORLIST(i,2);y=VECTORLIST(i,1);vx=VECTORLIST(i,4);vy=VECTORLIST(i,3);
% quiver(x,y,vx,vy,0.5,'LineWidth',vectorwidth,'Color','k','MaxHeadSize',vectorwidth);  
% set(gca,'XLim',[1 size(PLOTBIN,2)],'YLim',[1 size(PLOTBIN,1)],'YDir','reverse'); 
%  
%  %plot vectors of interest in separate figure
%  vectorfigure=figure;
%  hold on
%  quiver(0,0.5,N(1),N(2),0.5,'LineWidth',vectorwidth,'Color','r','MaxHeadSize',vectorwidth/2);
%  quiver(0,0.5,V(1),V(2),0.5,'LineWidth',vectorwidth,'Color','k','MaxHeadSize',vectorwidth/2);
%  set(gca,'XLim',[-1 1],'YLim',[-1 1],'YDir','reverse');   
%  