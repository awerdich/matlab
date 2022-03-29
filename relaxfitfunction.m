%function to fit ca uptake phase

function output=relaxfitfunction(P,X)
global Z1 V1
%fit parameter

output=exp(-(X-P(3))*(P(1)-V1)/Z1).*(P(1)/P(2)*sin(P(2)*(X-P(3)))+Z1*cos(P(2)*(X-P(3))))+P(4);
