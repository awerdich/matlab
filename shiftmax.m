DATACA=CA(:,1);
DATASA=-SA(:,1);
%maxima
%fit second order polynom for Ca
[val,idx]=max(DATACA);
DATAX=[idx-100:idx+1500]';
DATAY=DATACA(DATAX);
p=polyfit(DATAX,DATAY,2);
DATAFITY=polyval(p,DATAX);
[val,idxfit]=max(DATAFITY);idxfit=idxfit+DATAX(1)-1;
maxca=idxfit;
%fit for sarcomere length
[val,idx]=max(DATASA);
maxsa=idx;
%plot fit results
figure
plot(DATAX,DATAY);hold on
plot(DATAX,DATAFITY,'r');
plot(idxfit,polyval(p,idxfit),'r+')
figure
plot(DATACA,'k');hold on
plot(maxca,DATACA(maxca),'r+');
plot(DATASA*10,'b')
plot(maxsa,DATASA(maxsa)*10,'r+');
%figure delay
delayscans=maxsa-maxca;