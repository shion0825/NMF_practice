clear; close all; clc;

% kLen : ランクK
% en : 繰り返し回数N
kLen = 32;
en = 1024;
F = DGTtool(windowShift=512, windowLength=2048, FFTnum=2048, windowName="Hann");
NMF = NMF(kLen, en);

% xVec : 入力信号x
% xMat : 複素スペクトラムX
% xAbsMat : 振幅スペクトラムX
% xPhaseMat : 位相スペクトラムX
[xVec, fs] = audioread("in/t2in.wav");
xMat = F(xVec);
xAbsMat = abs(xMat);
xPhaseMat = angle(xMat);

% Euで計算
[yAbsMat, ~, ~, errVecEu] = NMF.calcNMF(xAbsMat, "Eu");
yMat = yAbsMat .* exp(1i * xPhaseMat);
yVec = F.pinv(yMat);
audiowrite("out/t2outEu.wav", yVec, fs);

% KLで計算
[yAbsMat, ~, ~, errVecKL] = NMF.calcNMF(xAbsMat, "KL");
yMat = yAbsMat .* exp(1i * xPhaseMat);
yVec = F.pinv(yMat);
audiowrite("out/t2outKL.wav", yVec, fs);

% ISで計算
[yAbsMat, ~, ~, errVecIS] = NMF.calcNMF(xAbsMat, "IS");
yMat = yAbsMat .* exp(1i * xPhaseMat);
yVec = F.pinv(yMat);
audiowrite("out/t2outIS.wav", yVec, fs);

semilogy(errVecEu);
hold on
semilogy(errVecKL);
semilogy(errVecIS);
legend("Eu", "KL", "IS");
hold off
