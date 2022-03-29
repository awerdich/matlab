%% open file
[avifile avipath]=uigetfile('*.avi','SELECT MOVIE');
mov=aviread([avipath,avifile]);
%% select angle
angle=input('