clc
clear
close all

basename = 'phs_test_';

N = fliplr([40, 80, 120, 160, 200, 240, 280, 320]);
M = 4*N;

Kmul = 8;

d = [3,5,7];
p = [3,5,7,9];
nm = [1.5, 2, 3];
g = [5, 10, 20];


leg = {'BSeuCall2Dbasket_FD', 'BSeuCall2Dbasket_RBFFDreg_phs', 'BSeuCall2Dbasket_RBFFDreg_smooth_phs', 'BSeuCall2Dbasket_RBFFDadap_phs'};

legN = numel(leg);
hN = numel(N);
