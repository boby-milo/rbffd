clc
clear
% close all

basename = 'phs_test_';

% N = fliplr([40, 80, 120, 160, 200, 240, 280, 320]);
% N = [40, 80, 120, 160, 200, 240, 280, 320, 360, 400];
% N = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,340,360,];
N = 40:20:400;
M = 10*N;

Kmul = 4;

d = 7;
p = 4;
nm = 3;
g = 5;


leg = {'BSeuCall2Dbasket_FD2', 'BSeuCall2Dbasket_RBFFDreg_phs', 'BSeuCall2Dbasket_RBFFDreg_smoothed_phs', 'BSeuCall2Dbasket_RBFFDadap_phs','BSeuCall2Dbasket_RBFFDadap_smoothed_phs','BSeuCall2Dbasket_RBFFDrepel_phs','BSeuCall2Dbasket_RBFFDrepel_smoothed_phs'};

legN = numel(leg);
hN = numel(N);
