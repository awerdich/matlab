%de-noising of the signal
%calculate default parameters
[C,L]=wavedec(A,3,'db10');

%default parameter for de-noising
[thr,sorh,keepapp]=ddencmp('den','wv',A);

%signal reconstruction
clean=wdencmp('gbl',C,L,'db10',3,thr,sorh,keepapp);