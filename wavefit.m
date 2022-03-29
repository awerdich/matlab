%function for fitting activation times
function output=wavefit(P,DATA)
X=DATA(:,1);
Y=DATA(:,2);
output=P(1)*X.^2+P(2)*Y.^2+P(3)*X.*Y+P(4)*X+P(5)*Y+P(6);