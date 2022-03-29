%function [CSCANINTERVALS]=correctseparation(WAVELENGTHS,SCANINTERVALS,PIXEL,PIXELINTERVAL)
%% correct SCANINTERVALS using fluorescence data if READY signal unavailable
CSCANINTERVALS=zeros(size(SCANINTERVALS));
i=1;XDATA=[];
while min(SCANINTERVALS(i,:))>0;%stop retrieving data when data for at least one wavelength is missing (e.g. end)
    for lambda=1:length(WAVELENGTHS)
        I=[SCANINTERVALS(i,2*lambda-1)+1;SCANINTERVALS(i,2*lambda)];
        %check for extreme fluorescence values
        X=(I(1):I(2))';Y=PIXEL(X);
        reference=PIXEL(round(mean(X)));
        fdev=15;%maximum deviation for fluorescence values
        
        %check if first value is out of line
        if abs(Y(1)-reference>fdev
            X(1)=[];Y(1)=[];
        else
            %if first value is ok, try previous one
            if Y(1)-1>0 && abs(PIXEL(Y(1)-1)-mean(Y))<fdev
                X=[X(1)-1;X];Y=PIXEL(X);
            end
        end

        if length(Y)>2
            %check last value if more than2 points left
            if abs(Y(end)-mean(Y))>fdev
                X(end)=[];Y(end)=[];
            end
        end
        %try next one
        k=X(end);
        while abs(PIXEL(k)-mean(Y))<fdev,k=k+1;end
        X=[X;[X(end):(k-1)]'];Y=PIXEL(X);
        CSCANINTERVALS(i,2*lambda-1)=X(1);
        CSCANINTERVALS(i,2*lambda)=X(end);
    end
    i=i+1;
end