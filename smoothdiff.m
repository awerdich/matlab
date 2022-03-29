%filter transients to calculate differential

function [D0,D1,D2]=smoothdiff(DATA,golaywindow)

%smoothing parameters
order=4; %smoothdiff: Order of polynomial fit
%golaywindow=31;%sent to Khaled on 10/6/2010 for cplexapd

%calculate filter coefficients
[b,g] = sgolay(order,golaywindow);   % Calculate S-G coefficients

%initialize results
HalfWin  = ((golaywindow+1)/2) -1;
D0=zeros(size(DATA)-HalfWin-1);
D1=zeros(size(DATA)-HalfWin-1);
D2=zeros(size(DATA)-HalfWin-1);

for i=1:size(DATA,2) 
    %extract transient
    T=DATA(:,i);

    for n = (golaywindow+1)/2:length(T)-(golaywindow+1)/2,
        % Zero-th derivative (smoothing only)
        SG0(n) =   dot(g(:,1), T(n - HalfWin: n + HalfWin));
  
        % 1st differential
        SG1(n) =   dot(g(:,2), T(n - HalfWin: n + HalfWin));
  
        % 2nd differential
        %SG2(n) = 2*dot(g(:,3)', MED_T(n - HalfWin: n + HalfWin))';
    end
    %Output results
    D0(:,i)=SG0';
    D1(:,i)=SG1';         
    %D2(:,i)=SG2';
end