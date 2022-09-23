%function to fit photobleaching by single exponential
%f(t)=A*exp(-t/tau)

function output=bleachfunction(P,X)

%fit parameter
A=P(1);
tau=P(2);

output=A*exp(-X/tau);