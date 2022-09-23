%three-point estimate of the derivative of columns in DATA
function OUTPUT=diff3(DATA)

OUTPUT=zeros(size(DATA));

for j=1:size(DATA,2)
    for i=2:size(DATA,1)-1
        OUTPUT(i,j)=(DATA(i+1,j)-DATA(i-1,j))/2;
    end
    %pad first and last point to retain original length of DATA
    OUTPUT(1,j)=OUTPUT(2,j);
    OUTPUT(end,j)=OUTPUT(end-1,j);
end

