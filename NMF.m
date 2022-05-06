function yMat = NMF(xMat, kLen, loop, div)
% [入力]
% xmat : 入力行列
% kLen : ランク
% loop : 更新回数
% div  : 距離の計算方法(Eu, KL, IS)

% [出力]
% yMat  : 分解後に再生成した行列

% wMat : 基底行列W
% hMat : 係数行列H
% yMat : 出力行列Y
% oMat : 全要素が1の行列O
[iLen, jLen] = size(xMat);
wMat = rand(iLen, kLen);
hMat = rand(kLen, jLen);
yMat = wMat * hMat;
oMat = ones(iLen, jLen);

% errVec : 誤差
errVec = zeros(loop, 1);
switch div
    case "Eu"
        for i = 1:loop
            wMat = wMat .* ((xMat * hMat') ./ (wMat * (hMat * hMat')));
            hMat = hMat .* ((wMat' * xMat) ./ ((wMat' * wMat) * hMat));
            yMat = wMat * hMat;
            errVec(i) = sum((xMat - yMat) .^ 2, "all");
        end
    case "KL"
        for i = 1:loop
            wMat = wMat .* (((xMat ./ yMat) * hMat') ./ (oMat * hMat'));
            yMat = wMat * hMat;
            hMat = hMat .* ((wMat' * (xMat ./ yMat)) ./ (wMat' * oMat));
            yMat = wMat * hMat;
            errVec(i) = sum(xMat .* log(xMat ./ yMat) - (xMat - yMat), "all");
        end
    case "IS"
        for i = 1:loop
            wMat = wMat .* (((xMat ./ (yMat .^ 2)) * hMat') ./ ((oMat ./ yMat) * hMat')) .^ 0.5;
            yMat = wMat * hMat;
            hMat = hMat .* ((wMat' * (xMat ./ (yMat .^ 2))) ./ (wMat' * (oMat ./ yMat))) .^ 0.5;
            yMat = wMat * hMat;
            errVec(i) = sum(xMat ./ yMat - log(xMat ./ yMat) - 1, "all");
        end
    otherwise
        disp("Skipped");
        return;
end

semilogy(errVec);
end