GL=A(:,2);PA=A(:,3);
%calculate offset
ogl=mean(GL(1:1000));opa=mean(PA(1:1000));
%subtract offsets
BGL=GL-ogl;BPA=PA-opa;
%determine normalization factors
fpa=700/max(BPA);
%normalize both transients with same factor
NGL=BGL*fpa;NPA=BPA*fpa;
%add offsets
CAPA=NPA+500;
CAGL=NGL+500*opa/og;