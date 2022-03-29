%return alphamatrix ALPHA for image I
function ALPHA=alphamatrix(I,level,transparency)

ALPHA=ones(size(I))*transparency;

for i=1:size(I,1)
    for j=1:size(I,2)
        if I(i,j)<level
            ALPHA(i,j)=0;
        end
    end
end