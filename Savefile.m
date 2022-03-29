%%Extract Values
N = zeros(51,1);
nu = input('heart #:');
name = DATABASE(nu).name;
N(1,:) = DATABASE(nu).age;
N(2,:) = DATABASE(nu).heart;
N(3,:) = DATABASE(nu).rate;
%N(4,:) = DATABASE(nu).av_maxrmse;
N(5,:) = DATABASE(nu).av_vel;
%N(6,:) = DATABASE(nu).av_velstd;
N(7,:) = DATABASE(nu).av_meanapd_ms20;
%N(8,:) = DATABASE(nu).av_meanvmax_s;
%N(9,:) = DATABASE(nu).av_meanvmaxsm_s;
%N(10,:) = DATABASE(nu).v_maxrmse;
N(11,:) = DATABASE(nu).v_vel;
%N(12,:) = DATABASE(nu).v_velstd;
N(13,:) = DATABASE(nu).v_meanapd_ms20;
%N(14,:) = DATABASE(nu).v_meanvmax_s;
%N(15,:) = DATABASE(nu).v_meanvmaxsm_s;
%N(16,:) = DATABASE(nu).a_maxrmse;
N(17,:) = DATABASE(nu).a_vel;
%N(18,:) = DATABASE(nu).a_velstd;
N(19,:) = DATABASE(nu).a_meanapd_ms20;
%N(20,:) = DATABASE(nu).a_meanvmax_s;
%N(21,:) = DATABASE(nu).a_meanvmaxsm_s;
%N(22,:) = DATABASE(nu).vic_maxrmse;
N(23,:) = DATABASE(nu).vic_vel;
%N(24,:) = DATABASE(nu).vic_velstd;
N(25,:) = DATABASE(nu).vic_meanapd_ms20;
%N(26,:) = DATABASE(nu).vic_meanvmax_s;
%N(27,:) = DATABASE(nu).vic_meanvmaxsm_s; 
%N(28,:) = DATABASE(nu).voc_maxrmse;
N(29,:) = DATABASE(nu).voc_vel;
%N(30,:) = DATABASE(nu).voc_velstd;
N(31,:) = DATABASE(nu).voc_meanapd_ms20;
%N(32,:) = DATABASE(nu).voc_meanvmax_s;
%N(33,:) = DATABASE(nu).voc_meanvmaxsm_s;
%N(34,:) = DATABASE(nu).aic_maxrmse;
N(35,:) = DATABASE(nu).aic_vel;
%N(36,:) = DATABASE(nu).aic_velstd;
N(37,:) = DATABASE(nu).aic_meanapd_ms20;
%N(38,:) = DATABASE(nu).aic_meanvmax_s;
%N(39,:) = DATABASE(nu).aic_meanvmaxsm_s;
%N(40,:) = DATABASE(nu).aoc_maxrmse;
N(41,:) = DATABASE(nu).aoc_vel;
%N(42,:) = DATABASE(nu).aoc_velstd;
N(43,:) = DATABASE(nu).aoc_meanapd_ms20;
%N(44,:) = DATABASE(nu).aoc_meanvmax_s;
%N(45,:) = DATABASE(nu).aoc_meanvmaxsm_s;
%N(46,:) = DATABASE(nu).o_maxrmse;
%N(47,:) = DATABASE(nu).o_vel;
%N(48,:) = DATABASE(nu).o_velstd;
%N(49,:) = DATABASE(nu).o_meanapd_ms20;
%N(50,:) = DATABASE(nu).o_meanvmax_s;
%N(51,:) = DATABASE(nu).o_meanvmaxsm_s;
%%
%%exporting and saving csv file
%if exist('fn') == 1; 
%    dlmwrite(fn,N,'-append');
%else    
    
fn = [name, '.csv'];
dlmwrite(fn,N),',';
close all;
%end