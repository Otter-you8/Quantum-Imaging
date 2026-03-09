%% 各画素における強度相関関数を用いて、minus-coordinateを求める
clear; clc; close all;
fprintf('"Minus-coordinate"\n');
fprintf('-----------------------------------------------------\n');

%% 設定項目
% 【追加】GUIでフォルダのパスを指定
folderPath = uigetdir(pwd, 'TIFFファイルが含まれるフォルダを選択してください');
if folderPath == 0
    fprintf('フォルダの選択がキャンセルされました。処理を終了します。\n');
    return; % キャンセル時は処理を終了
end

% フォルダ内のtifファイルを取得
fileList = dir(fullfile(folderPath, '*.tif'));
dataSet = length(fileList);

if dataSet == 0
    fprintf('選択したフォルダにTIFFファイルが見つかりませんでした。処理を終了します。\n');
    return;
end

% 画像サイズ設定
fprintf('最初のファイルから画像サイズを取得しています...\n');
firstFilePath = fullfile(folderPath, fileList(1).name);
temp_array = readFile(firstFilePath);
[height, width, ~] = size(temp_array);
num_pixels = height * width;
clear temp_array; % メモリ節約のため一時変数をクリア
fprintf('画像サイズ: Width = %d, Height = %d\n', width, height);

% その他の設定
binary = false; % ※可能であれば後日true（閾値処理）での検証も推奨します
threshold_value = 540;
threshold_max = 20000;

max_total_frames = 1000000;  
margin = 10;

%% 相関計算
sum_corr1 = zeros(num_pixels, num_pixels, 'double'); 
sum_corr2 = zeros(num_pixels, num_pixels, 'double'); 

total_pairs_corr1 = 0; % 第1項の計算に使われた有効ペア数
total_pairs_corr2 = 0; % 第2項の計算に使われた有効ペア数

all_sum = zeros(height, width);
current_total_frames = 0;

tStart = tic;
for dataSetNumber = 1:dataSet
    if current_total_frames >= max_total_frames
        fprintf('指定枚数 (%d) に達したため読み込みを終了します。\n', max_total_frames);
        break; 
    end

    fprintf('Dataset : %d / %d\n', dataSetNumber, dataSet);
    filePath = fullfile(folderPath, fileList(dataSetNumber).name);
    
    fprintf('  Loading files...\n');
    original_array = readFile(filePath);
    frames_in_file = size(original_array, 3);

    if current_total_frames + frames_in_file > max_total_frames
        needed_frames = max_total_frames - current_total_frames;
        if needed_frames > 0
            original_array = original_array(:, :, 1:needed_frames);
            fprintf('  (Limit reached: Using first %d frames of this file)\n', needed_frames);
        else
            break;
        end
    end

    if size(original_array, 3) < 2
        continue;
    end

    % 宇宙線の影響があるフレームを除去
    max_intensity = squeeze(max(original_array, [], [1,2]));
    is_valid_frame = max_intensity < threshold_max;
    original_array = original_array(:, :, is_valid_frame);
    
    removed_frame = find(~is_valid_frame);
    if ~isempty(removed_frame)
        fprintf('  除去されたフレーム数: %d\n', length(removed_frame));
    end

    if binary == true
        original_array = single(original_array > threshold_value);
    end

    N = size(original_array, 3);
    if N < 2
        fprintf('  Skipping: Not enough frames after filtering (N=%d).\n', N);
        continue; 
    end
    
    current_total_frames = current_total_frames + N;
    all_sum = all_sum + sum(original_array, 3);

    fprintf('  Calculation of intensity correlation...\n');

    % 1D配列化 (hw x N)
    reshaped_array = reshape(permute(original_array,[2,1,3]), num_pixels, N); 
    
    % 第１項（同じフレーム間の積の総和）
    fprintf('   corr1.\n');
    sum_corr1 = sum_corr1 + (reshaped_array * reshaped_array');
    total_pairs_corr1 = total_pairs_corr1 + N;

    % 第２項（連続するフレーム間の積の総和）
    fprintf('   corr2.\n');
    reshaped_array_1 = reshaped_array(:, 1:N-1);
    reshaped_array_2 = reshaped_array(:, 2:N);
    sum_corr2 = sum_corr2 + (reshaped_array_1 * reshaped_array_2' + reshaped_array_2 * reshaped_array_1');
    total_pairs_corr2 = total_pairs_corr2 + 2*(N-1);
end
CalculationTime = toc(tStart);
fprintf('Total calculation time: %.2f seconds\n', CalculationTime);

%% 全フレームでの平均化と減算
intensity_corr1 = sum_corr1 / total_pairs_corr1;
intensity_corr2 = sum_corr2 / total_pairs_corr2;

% 第1項から第2項（偶発同時計数）を減算
intensityCorr_all = intensity_corr1 - intensity_corr2;

% 相関マップの形になおす [height, width, hw]
corrMap_all = reshape(intensityCorr_all', [width, height, num_pixels]);
corrMap_all = permute(corrMap_all, [2, 1, 3]);

%% 自己相関成分（対角成分）の補間
for y = 1:height
    for x = 1:width
        idx = (y-1)*height + x;
        if x > 1 && x < width
            corrMap_all(y, x, idx) = (corrMap_all(y, x-1, idx) + corrMap_all(y, x+1, idx)) / 2;
        end
    end
end

%% クロストークの補正
delta = 3;
padded_intensityCorr = padarray(corrMap_all, [delta, delta, 0], 0, 'both');
corrected_intensityCorr = padded_intensityCorr; % 書き込み用の配列を分離

index = 1;
for y_value = delta+1 : delta+height
    for x_value = delta+1 : delta+width
        temp_old = padded_intensityCorr(:,:,index); % 読み込み専用
        temp_new = temp_old; % 書き込み用
        
        for dp = -3 : 3
            if (y_value + dp >= 1 && y_value + dp <= height+2*delta && x_value - 1 >= 1 && x_value + 1 <= width+2*delta)
                temp_new(y_value, x_value + dp) = (temp_old(y_value-1, x_value + dp) + temp_old(y_value+1, x_value + dp)) / 2;
            end
        end
        temp_new(:, x_value) = (temp_old(:, x_value-1) + temp_old(:, x_value+1)) / 2;
        
        corrected_intensityCorr(:,:,index) = temp_new;
        index = index + 1;
    end
end

% 0パディング部分を除去
intensityCorr_all = corrected_intensityCorr(delta+1:delta+height, delta+1:delta+width, :);

%% minus-coordinateへの投影
minusCoordinate = zeros(2*height, 2*width);

for x = 1:width
    for y = 1:height
        idx = (y-1)*height + x;
        currentMap = intensityCorr_all(:,:,idx);

        y_start = height + 1 - (y-1);
        y_end = y_start + (height-1);
        x_start = width + 1 - (x-1);
        x_end = x_start + (width-1);

        minusCoordinate(y_start:y_end, x_start:x_end) = minusCoordinate(y_start:y_end, x_start:x_end) + currentMap;
    end
end

[h, w] = size(minusCoordinate);
cx = ceil((w+1)/2);
cy = ceil((h+1)/2);

% smearingの補正（投影図の中心）
minusCoordinate(cy,:) = (minusCoordinate(cy-1,:) + minusCoordinate(cy+1,:)) / 2;

% 切り取り
dy = floor(height/2);
dx = floor(width/2);
minusCoordinate_2 = minusCoordinate(cy-dy:cy+dy, cx-dx:cx+dx);

%% 結果の表示
figure(1);
imagesc(minusCoordinate_2);
axis equal tight
set(gcf, 'Position', [500, 400, 600, 300]);
title('Minus coordinate');
xlabel('x1-x2');
ylabel('y1-y2');
colorbar;
colormap("jet");

figure(2);
imagesc(all_sum);
axis equal tight
title('All Sum Image');
colormap("gray"); colorbar;

%% ファイルを読み込む関数
function spdc_photons0 = readFile(tiff_photonsFile)
  spdc_photons0 = tiffreadVolume(tiff_photonsFile);
  spdc_photons0 = single(spdc_photons0);
end