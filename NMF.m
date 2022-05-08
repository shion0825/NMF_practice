classdef NMF
    properties
        kLen {mustBeInteger, mustBePositive}
        loop {mustBeInteger, mustBePositive}
    end

    methods (Access = public)
        function obj = NMF(kLen, loop)
            % [入力]
            % kLen : 分解ランク
            % loop : 更新回数
            % [出力]
            % obj : NMFクラスのインスタンス
            obj.kLen = kLen;
            obj.loop = loop;
        end
    end

    methods (Access = public)
        function [yMat, wMat, hMat, errVec] = calcNMF(obj, xMat, div)
            % [入力]
            % xMat : 入力行列X
            % div : 距離の計算方法(Eu, KL, IS)
            % [出力]
            % yMat : 出力行列Y
            % wMat : 基底行列W
            % hMat : 係数行列H
            % errVec : 誤差
            [iLen, jLen] = size(xMat);
            wMat = rand(iLen, obj.kLen);
            hMat = rand(obj.kLen, jLen);
            yMat = wMat * hMat;
            oMat = ones(iLen, jLen);
            errVec = zeros(obj.loop, 1);

            switch div
                case "Eu"
                    for i = 1:obj.loop
                        wMat = wMat .* ((xMat * hMat') ./ (yMat * hMat'));
                        yMat = wMat * hMat;
                        hMat = hMat .* ((wMat' * xMat) ./ (wMat' * yMat));
                        yMat = wMat * hMat;
                        errVec(i) = sum((xMat - yMat) .^ 2, "all");
                    end
                case "KL"
                    for i = 1:obj.loop
                        wMat = wMat .* (((xMat ./ yMat) * hMat') ./ (oMat * hMat'));
                        yMat = wMat * hMat;
                        hMat = hMat .* ((wMat' * (xMat ./ yMat)) ./ (wMat' * oMat));
                        yMat = wMat * hMat;
                        errVec(i) = sum(xMat .* log(xMat ./ yMat) - (xMat - yMat), "all");
                    end
                case "IS"
                    for i = 1:obj.loop
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
        end

        function [yMat, gMat1, gMat2, errVec] = calcActivationMat2(obj, mMat, wMat1, wMat2, div)
            % [入力]
            % mMat : 混合行列X
            % wMat1 : 基底行列W1
            % wMat2 : 基底行列W2
            % div : 距離の計算方法(Eu, KL, IS)
            % [出力]
            % yMat : 出力行列Y
            % gMat1 : 係数行列G1
            % gMat2 : 係数行列G2
            % errVec : 誤差
            [iLen, jLen] = size(mMat);
            gMat1 = rand(obj.kLen, jLen);
            gMat2 = rand(obj.kLen, jLen);
            yMat = wMat1 * gMat1 + wMat2 * gMat2;
            oMat = ones(iLen, jLen);
            errVec = zeros(obj.loop, 1);

            switch div
                case "Eu"
                    for i = 1:obj.loop
                        gMat1 = gMat1 .* ((wMat1' * mMat) ./ (wMat1' * yMat));
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        gMat2 = gMat2 .* ((wMat2' * mMat) ./ (wMat2' * yMat));
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        errVec(i) = sum((mMat - yMat) .^ 2, "all");
                    end
                case "KL"
                    for i = 1:obj.loop
                        gMat1 = gMat1 .* ((wMat1' * (mMat ./ yMat)) ./ (wMat1' * oMat));
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        gMat2 = gMat2 .* ((wMat2' * (mMat ./ yMat)) ./ (wMat2' * oMat));
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        errVec(i) = sum(mMat .* log(mMat ./ yMat) - (mMat - yMat), "all");
                    end
                case "IS"
                    for i = 1:obj.loop
                        gMat1 = gMat1 .* ((wMat1' * (mMat ./ (yMat .^ 2))) ./ (wMat1' * (oMat ./ yMat))) .^ 0.5;
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        gMat2 = gMat2 .* ((wMat2' * (mMat ./ (yMat .^ 2))) ./ (wMat2' * (oMat ./ yMat))) .^ 0.5;
                        yMat = wMat1 * gMat1 + wMat2 * gMat2;
                        errVec(i) = sum(mMat ./ yMat - log(mMat ./ yMat) - 1, "all");
                    end
                otherwise
                    disp("Skipped");
                    return;
            end
        end
    end
end