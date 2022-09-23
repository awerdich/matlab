%Extract Variables
N = zeros(51,1);
nu = input('heart #:');
name = FDATABASE(nu).name;
%N(1,:) = FDATABASE(nu).age;
%N(2,:) = FDATABASE(nu).heart;
%N(3,:) = FDATABASE(nu).rate;
N(4,:) = FDATABASE(nu).av_diast;
N(5,:) = FDATABASE(nu).av_amp;
N(6,:) = FDATABASE(nu).av_durms50;
N(7,:) = FDATABASE(nu).v_diast;
N(8,:) = FDATABASE(nu).v_amp;
N(9,:) = FDATABASE(nu).v_durms50;
N(10,:) = FDATABASE(nu).a_diast;
N(11,:) = FDATABASE(nu).a_amp;
N(12,:) = FDATABASE(nu).a_durms50;
%N(13,:) = FDATABASE(nu).vic_diast;
%N(14,:) = FDATABASE(nu).vic_amp;
%N(15,:) = FDATABASE(nu).vic_durms50;
%N(16,:) = FDATABASE(nu).voc_diast;
%N(17,:) = FDATABASE(nu).voc_amp;
%N(18,:) = FDATABASE(nu).voc_durms50;
%N(19,:) = FDATABASE(nu).aic_diast;
%N(20,:) = FDATABASE(nu).aic_amp;
%N(21,:) = FDATABASE(nu).aic_durms50;
%N(22,:) = FDATABASE(nu).aoc_diast;
%N(23,:) = FDATABASE(nu).aoc_amp;
%N(24,:) = FDATABASE(nu).aoc_durms50;
%N(25,:) = FDATABASE(nu).o_diast;
%N(26,:) = FDATABASE(nu).o_amp;
%N(27,:) = FDATABASE(nu).o_durms50;
%%
%%exporting and saving csv file
%if exist('fn') == 1; 
%    dlmwrite(fn,N,'-append');
%else    
    
fn = [name, '.csv'];
dlmwrite(fn,N),',';
close all;
%end