%function to fit ca uptake phase

function output=uptakefitfunction(P,X)

output=abs(P(1)*exp(-X/P(2)));