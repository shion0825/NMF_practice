clear; close all; clc;

% kLen : ランクK
% en : 繰り返し回数N
kLen = 16;
en = 256;
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
yMat = NMF(xAbsMat, kLen, en, "KL");
yVec = F.pinv(yMat);
%yVec = 0.8 * yVec / max(yVec, [], "all");
audiowrite("out/out.wav", yVec, fs);