%% FFTを用いて、位置相関を求める

clear; clc; close all;
fprintf('"Position Correlatin Analysis using FFT"\n');
fprintf('-----------------------------------------------------\n');

%% パラメータ設定
% フォルダのパス設定
folderPath = 'C:\Users\cs13\Pictures\EMCCD\Correlation\260123_20ms_NoBFP';
%フォルダ内のtifファイルを取得
fileList = dir(fullfile(folderPath, '*.tif'));
% 計算するセット数
dataSet = length(fileList);

% 上限閾値処理（異常イベントの除去）
threshold_max = 30000;

% 計算フレーム数の制限
max_total_frames = 1000000;  % ここで指定（例：1000枚で計算を打ち切る）
current_total_frames = 0; % 現在の積算枚数カウンタ

% プロット設定
cropSize = 50; % 相関マップとして切り出すサイズ

% 相関計算に使ったフレーム数
valid_frames_count = 0;

% 移動平均のフレーム数
half_window = 1000;

%% 相関計算

if dataSet > 0
    firstFile = fullfile(folderPath, fileList(1).name);
    info = imfinfo(firstFile);
    Height = info(1).Height;
    Width = info(1).Width;
else
    error('フォルダにtifファイルが見つかりません。');
end

num_frames = zeros(1,dataSet); % 枚数を格納する配列

padH = 2*Height;
padW = 2*Width;

totalAutoCorr = zeros(padH, padW);
totalCrossCorr = zeros(padH, padW);

All_Sum = zeros(Height, Width);

tStart = tic;
for dataSetNumber = 1:dataSet
    
    if current_total_frames >= max_total_frames
        fprintf('指定枚数 (%d) に達したため読み込みを終了します。\n', max_total_frames);
        break; 
    end

    fprintf('Dataset : %d\n', dataSetNumber);
    filePath = fullfile(folderPath, fileList(dataSetNumber).name);
   
    % tifファイル読み込み
    fprintf(' Loading files...\n');
    clear original_array
    original_array = readFile(filePath);

    [h,w,frames_in_file] = size(original_array);
    if h ~= Height || w ~= Width
        fprintf(' Error: Image size mismatch. Skipping file. \n');
        continue;
    end

    if current_total_frames + frames_in_file > max_total_frames
        needed_frames = max_total_frames - current_total_frames;
        if needed_frames > 0
            original_array = original_array(:, :, 1:needed_frames);
            fprintf('  (Limit reached: Using first %d frames of this file)\n', needed_frames);
        else
            % 既に満たしている場合はスキップして終了
            break;
        end
    end

    % フレーム数が計算（相関計算には最低2枚必要）に足りない場合の安全策
    if frames_in_file < 2
        fprintf('  Skipping: Not enough frames (%d) for correlation.\n', frames_in_file);
        continue;
    end
    
    % 総フレーム数
    current_total_frames = current_total_frames + frames_in_file;
    num_frames(dataSetNumber) = current_total_frames;

    % 各フレームの最大値に対して閾値以下のフレームを特定する（宇宙線の影響があるフレームを除去）
    clear max_intensity is_valid_frame removed_frame
    max_intensity = squeeze(max(original_array, [], [1,2]));
    is_valid_frame = max_intensity < threshold_max;
    original_array = original_array(:, :, is_valid_frame);
    
    removed_frame = find(~is_valid_frame);

    if ~isempty(removed_frame)
        fprintf(' 除去されたフレーム: %d\n', removed_frame);
    end

    N = size(original_array, 3);

    if N < 2
        fprintf('  Skipping: Not enough frames after filtering (N=%d).\n', N);
        continue; % このデータセット(ファイル)の処理を飛ばして次へ行く
    end

    All_Sum = All_Sum + sum(original_array, 3);
    
    % 平均減算
%     mean_Image = mean(original_array,3);
%     fluctuations = original_array - mean_Image;

    background_stack = movmean(original_array, [half_window, half_window], 3);
    fluctuations = original_array - background_stack;

%     fluctuations = original_array;
    
    % FFTを用いた相関計算
    sumAutoCorr = zeros(padH, padW); % 同じフレームでの自己相関
    sumCrossCorr = zeros(padH, padW); % フレームをずらした交互相関

    % プログレスバー（ファイルごと）
    hWait = waitbar(0, sprintf('Dataset %d Calculation...', dataSetNumber));

    for n = 1:N-2
        frame1 = fluctuations(:,:,n);
        frame2 = fluctuations(:,:,n+2);

        frame1 = frame1 - mean(frame1, 1);
        frame2 = frame2 - mean(frame2, 1);
        
        F1 = fft2(frame1, padH, padW);
        F2 = fft2(frame2, padH, padW);

        % A. 自己相関
        Corr1 = ifft2(abs(F1).^2);
        sumAutoCorr = sumAutoCorr + real(Corr1);
        
        % B. 1フレームずらした交互相関
        Corr2 = ifft2(F1 .* conj(F2));
        sumCrossCorr = sumCrossCorr + real(Corr2);

        if mod(n, 100) == 0
            waitbar(n / (N-2), hWait);
        end
    end
    close(hWait);

    totalAutoCorr = totalAutoCorr + sumAutoCorr;
    totalCrossCorr = totalCrossCorr + sumCrossCorr;
    valid_frames_count = valid_frames_count + (N-2);


end
CalculationTime = toc(tStart);

avgAutoCorr = totalAutoCorr / valid_frames_count;
avgCrossCorr = totalCrossCorr / valid_frames_count;

avgAutoCorr = fftshift(avgAutoCorr);
avgCrossCorr = fftshift(avgCrossCorr);

finalCorr = avgAutoCorr - avgCrossCorr;

% 表示範囲のためにクロップ (中心部分のみ切り出し)
cX = padW / 2 + 1;
cY = padH / 2 + 1;
% finalCorr(cY, cX) = NaN;
% finalCorr(cY, cX - 1) = NaN;
% finalCorr(cY, cX + 1) = NaN;

finalCorr(cY,:) = (finalCorr(cY+1,:) + finalCorr(cY-1,:)) / 2;

rangeX = cX - cropSize : cX + cropSize;
rangeY = cY - cropSize : cY + cropSize;

finalCorrCropped = finalCorr(rangeY, rangeX);
autoCorrCropped = avgAutoCorr(rangeY, rangeX);
bgCorrCropped = avgCrossCorr(rangeY, rangeX);

%% 結果のプロット

figure(1);
imagesc(autoCorrCropped);
title('Raw Autocorrelation (Signal + Noise)');
axis image;  colorbar;

figure(2);
imagesc(bgCorrCropped);
title('Background Correlation (Frame n & n+1)');
axis image;  colorbar;

figure(3);
imagesc(finalCorrCropped);
title('Subtracted Result (Quantum Correlation)');
axis image;  colorbar;

figure(4);
imagesc(All_Sum);
title('All Sum Image');
axis image;  colormap("gray");  colorbar;

% ピークの確認
disp('解析完了。Figureを確認してください。');


%% ファイルを読み込む関数
function spdc_photons0 = readFile(tiff_photonsFile)
  spdc_photons0 = tiffreadVolume(tiff_photonsFile);
  spdc_photons0 = single(spdc_photons0);
  % spdc_photons0 = spdc_photons0 / 255;
end
