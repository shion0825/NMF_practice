clear; close all; clc;

% kLen : ランクK
% en : 繰り返し回数N
kLen = 2;
en = 128;
F = DGTtool(windowShift=512, windowLength=2048, FFTnum=2048, windowName="Hann");
NMF = NMF(kLen, en);

xMat = [1, 2, 3; 3, 6, 9; 4, 5, 6; 6, 7.5, 9]';
[~, ~, ~, errVecEu] = NMF.calcNMF(xMat, "Eu");
[~, ~, ~, errVecKL] = NMF.calcNMF(xMat, "KL");
[~, ~, ~, errVecIS] = NMF.calcNMF(xMat, "IS");

semilogy(errVecEu);
hold on
semilogy(errVecKL);
semilogy(errVecIS);
legend("Eu", "KL", "IS");
hold off