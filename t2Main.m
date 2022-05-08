clear; close all; clc;

% kLen : ランクK
% en : 繰り返し回数N
kLen = 32;
en = 1024;
F = DGTtool(windowShift=512, windowLength=2048, FFTnum=2048, windowName="Hann");

% xVec : 入力信号x
% xMat : 複素スペクトラムX
% xAbsMat : 振幅スペクトラムX
% xPhaseMat : 位相スペクトラムX
[xVec, fs] = audioread("in/in.wav");
xMat = F(xVec);
xAbsMat = abs(xMat);
xPhaseMat = angle(xMat);

% KLで計算(Eu, ISも可)
yAbsMat = NMF(xAbsMat, kLen, en, "Eu");
yMat = yAbsMat .* exp(1i * xPhaseMat);
yVec = F.pinv(yMat);
audiowrite("out/out.wav", yVec, fs);