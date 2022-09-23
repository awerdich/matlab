beforemax=1000;aftermax=1000;threshold=0.95;scanrate=2000;
TA=T;
BASELINE=[1:100];
%re-normalize data before alignment
TA=TA-mean(TA(BASELINE));
TA=TA/max(TA);

%determine index when threshold is reached from max
[val,idx]=max(TA);
j=idx;while(TA(j))>threshold;j=j-1;end;
A=TA(j-beforemax:j+aftermax-1);

%correct offset and normalize
%filter upstrokes
DA=diff3(A)*scanrate;
FA=filterdata(scanrate,60,DA);%low pass filter differential
DA=FA/max(FA);
