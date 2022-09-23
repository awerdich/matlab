%calculate divergence of vectorfield (VX,VY)
%VX(x,y): X-component of vectorfield
%VY(x,y): Y-component of vectorfield

function DN = cvdivergence(VX,VY)
%% run function as script
%define velocity vectors
VNX = zeros(size(VX));
VNY = zeros(size(VY));

%normalize vectors
for i=1:size(VX,1)
    for j=1:size(VX,2)
        V = [VX(i,j),VY(i,j)];
        vel = norm(V,2);
        %cannot divide by 0
        if vel>0
            VN = V/vel;
        else
            VN = [0,0];
        end
        VNX(i,j)=VN(1);
        VNY(i,j)=VN(2);
    end
end
      
%define coordinates of the velocity field
X=[1:size(VX,2)]';
Y=[1:size(VX,1)]';
[XM,YM]=meshgrid(X,Y);

%overall divergence
D=divergence(XM,YM,VNX,VNY);

%normalize D in [-1, +1]
minD = min(D(:));
maxD = max(D(:));
DN = zeros(size(D));
R=[-1,+1];
for i=1:size(D,1)
    for j=1:size(D,2)
        if abs(VX(i,j)) > 0
            %normalize to [0,1]
            DZ = (D(i,j) - minD)/(maxD-minD);
            %scale into R
            DN(i,j) = (DZ * (R(2)-R(1))) + R(1);
        end
    end
end
