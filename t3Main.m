clear; close all; clc;

% パス追加
addpath("./bss_eval");

% kLen : ランクK
% en : 繰り返し回数N
kLen = 64;
en = 512;
F = DGTtool(windowShift=512, windowLength=2048, FFTnum=2048, windowName="Hann");
NMF = NMF(kLen, en);

% xVec : 学習信号x
% sVec : 入力信号x
% mVec : 混合信号m
% xMat : 複素スペクトラムX
% mMat : 複素スペクトラムM
% xAbsMat : 振幅スペクトラムX
% mAbsMat : 振幅スペクトラムM
% mPhaseMat : 位相スペクトラムM
[xVec1, ~] = audioread("in/t3base1.wav");
[xVec2, ~] = audioread("in/t3base2.wav");
[sVec1, ~] = audioread("in/t3in1.wav");
[sVec2, fs] = audioread("in/t3in2.wav");
mVec = sVec1 + sVec2;
xMat1 = F(xVec1);
xMat2 = F(xVec2);
mMat = F(mVec);
xAbsMat1 = abs(xMat1);
xAbsMat2 = abs(xMat2);
mAbsMat = abs(mMat);
mPhaseMat = angle(mMat);

% Wを学習
[~, wMat1, ~, ~] = NMF.calcNMF(xAbsMat1, "KL");
[~, wMat2, ~, ~] = NMF.calcNMF(xAbsMat2, "KL");

% Gを計算
[yAbsMat, gMat1, gMat2, ~] = NMF.calcActivationMat2(mAbsMat, wMat1, wMat2, "KL");

% Yを出力
yMat = yAbsMat .* exp(1i * mPhaseMat);
yVec = F.pinv(yMat);
audiowrite("out/t3mix.wav", yVec, fs);

% W*Gを出力
wgAbsMat1 = wMat1 * gMat1;
wgAbsMat2 = wMat2 * gMat2;
wgMat1 = ((wgAbsMat1 .^ 2) ./ (wgAbsMat1 .^ 2 + wgAbsMat2 .^ 2)) .* mMat;
wgMat2 = ((wgAbsMat2 .^ 2) ./ (wgAbsMat1 .^ 2 + wgAbsMat2 .^ 2)) .* mMat;
wgVec1 = F.pinv(wgMat1);
wgVec2 = F.pinv(wgMat2);
audiowrite("out/t3sep1.wav", wgVec1, fs);
audiowrite("out/t3sep2.wav", wgVec2, fs);

% SDR計算
vLen = size(sVec1, 1);
[SDR, ~, ~] = bss_eval_sources([wgVec1(1:vLen), wgVec2(1:vLen)]', [sVec1, sVec2]');
disp(SDR);
