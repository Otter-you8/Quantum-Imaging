%% 各画素における強度相関関数を用いて、minus-coordinateを求める
clear; clc; close all;
fprintf('"Minus-coordinate"\n');
fprintf('-----------------------------------------------------\n');
% addpath('C:\MATLAB\intensity correlation\Functions\');

%% 設定項目

% スクリプトを実行して設定変数をロード
params;

%% 相関計算

% 強度相関関数
intensity_corr1 = zeros(height^2,width^2); % 第1項
intensity_corr2 = zeros(height^2,width^2); % 第2項

num_frames = zeros(1,dataSet); % 枚数を格納する配列

all_sum = zeros(height,width);
numFrames = 0;
figNumber = 1;

tStart = tic;
for dataSetNumber = 1:dataSet
    
    if current_total_frames >= max_total_frames
        fprintf('指定枚数 (%d) に達したため読み込みを終了します。\n', max_total_frames);
        break; 
    end

    fprintf('Dataset : %d / %d \n', dataSetNumber,dataSet);
    filePath = fullfile(folderPath, fileList(dataSetNumber).name);
   
    % tifファイル読み込み
    fprintf(' Loading files...\n');
    clear original_array
    original_array = readFile(filePath);

    frames_in_file = size(original_array, 3);

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

    N = size(original_array,3);

    % 総フレーム数
    current_total_frames = current_total_frames + N;
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

    if binary == true
        original_array = single(original_array > threshold_value);
    end

    N = size(original_array,3);

    if N < 2
        fprintf('  Skipping: Not enough frames after filtering (N=%d).\n', N);
        continue; % このデータセット(ファイル)の処理を飛ばして次へ行く
    end
    
    %総和画像
    all_sum = all_sum + sum(original_array,3);

    fprintf(' Calculation of intensity correlation...\n');

    % 第１項
    fprintf('  corr1.\n');
    
    reshaped_array = reshape(permute(original_array,[2,1,3]), 1,height*width,[]); 
    reshaped_array = squeeze(reshaped_array);
    intensity_corr1 = intensity_corr1 + (reshaped_array * reshaped_array')/N;
    intensity_corr1(1:size(intensity_corr1,1)+1:end) = 0; % 自己相関は０
    

    % 第２項
    fprintf('  corr2.\n');
    
    reshaped_array_1 = original_array(:,:,1:N-1);
    reshaped_array_1 = reshape(permute(reshaped_array_1,[2,1,3]), 1,height*width,[]);
    reshaped_array_1 = squeeze(reshaped_array_1);

    reshaped_array_2 = original_array(:,:,2:N);
    reshaped_array_2 = reshape(permute(reshaped_array_2,[2,1,3]), 1,height*width,[]);
    reshaped_array_2 = squeeze(reshaped_array_2);

    intensity_corr2 = intensity_corr2 + (reshaped_array_1*reshaped_array_2' + reshaped_array_2*reshaped_array_1')/(2*(N-1));
    intensity_corr2(1:size(intensity_corr2,1)+1:end) = 0;
    

end
CalculationTime = toc(tStart);

% 最終計算結果
intensityCorr_all = intensity_corr1 - intensity_corr2; % ここでは(height×width)×(height×width) 配列
intensityCorr_all(1:size(intensityCorr_all,1)+1:end) = 0;

% 相関マップの形になおす
corrMap_all = reshape(intensityCorr_all',[width,height,height*width]);
corrMap_all = permute(corrMap_all,[2,1,3]);

% 補正前の強度相関関数
%intensityCorr_all = intensity_corr1/(N*dataSet) - intensity_corr2/(N*dataSet-1);
intensityCorr = corrMap_all;
%% 

% すべてのページに０パディング
delta = 3;
padded_intensityCorr = padarray(intensityCorr, [delta, delta, 0], 0, 'both');

% クロストークの補正
index = 1;
for y_value = delta+1 : delta+height
    for x_value = delta+1 : delta+width
        temp = padded_intensityCorr(:,:,index);
        % クロストーク補正
        
        for dp = -3 : 3
            % y_value + dpが1以上height以下、x_value - 1が1以上、x_value + 1がwidth以下であることを確認
            if (y_value + dp >= 1 && y_value + dp <= height+2*delta && x_value - 1 >= 1 && x_value + 1 <= width+2*delta)
                temp(y_value, x_value + dp) = (temp(y_value-1, x_value + dp) + temp(y_value+1, x_value + dp)) / 2;
            end
        end
        
       temp(:,x_value) = (temp(:,x_value-1) + temp(:,x_value+1))/2;
       padded_intensityCorr(:,:,index) = temp;
       index = index + 1;
    end
end

% 0パディング部分を除去
intensityCorr_all = padded_intensityCorr(delta+1:delta+height,delta+1:delta+width,:);


% minus-coordinate用の配列
minusCoordinate = zeros(2*height, 2*width);
num_pixels = height*width;

for x = 1:width
    for y = 1:height
        % 相関マップが格納されているページを求める
        idx = (y-1)*height + x;
        currentMap = intensityCorr_all(:,:,idx);

        % fprintf('progress >>> (%d, %d)', x,y);
        % minus-coordinate用の配列におけるスタート位置を求める
        y_start = height + 1 - (y-1);
        y_end = y_start + (height-1);
        x_start = width + 1 - (x-1);
        x_end = x_start + (width-1);
        % fprintf('    y : %d ~ %d, x : %d ~ %d\n', y_start,y_end, x_start,x_end);

        % minus-coordinate配列を更新
        minusCoordinate(y_start:y_end, x_start:x_end) = minusCoordinate(y_start:y_end, x_start:x_end) + currentMap;

    end
end

%{
for idx = 1:num_pixels
    currentMap = intensityCorr_all(:,:,idx);

    % 相関マップのインデックス番号から、どの画素との相関マップかを求める
    y = ceil(idx / height);      % y座標
    x = mod(idx - 1, width) + 1; % x座標
    fprintf('progress >>> (%d, %d)', x,y);

    % minus-coordinate用の配列におけるスタート位置を求める
    y_start = height + 1 - (y-1);
    y_end = y_start + (height-1);
    x_start = width + 1 - (x-1);
    x_end = x_start + (width-1);
    fprintf('    y : %d ~ %d, x : %d ~ %d\n', y_start,y_end, x_start,x_end);

    % minus-coordinate配列を更新
    minusCoordinate(y_start:y_end, x_start:x_end) = minusCoordinate(y_start:y_end, x_start:x_end) + currentMap;

end
%}

[h, w] = size(minusCoordinate);
cx = ceil((w+1)/2);
cy = ceil((h+1)/2);

% smearingの補正
minusCoordinate(cy,:) =  (minusCoordinate(cy-1,:) + minusCoordinate(cy+1,:))/2;


% 切り取り
% p2 = size(minusCoordinate,1)/4;
% minusCoordinate_2 = minusCoordinate(cy-height/2:cy+height/2, cx-width/2:cx+width/2); % (heigh+1)x(width+1)

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

figure(2);
imagesc(all_sum);
axis equal tight
title('All Sum Image');
colormap("gray");  colorbar;


%% ファイルを読み込む関数
function spdc_photons0 = readFile(tiff_photonsFile)
  spdc_photons0 = tiffreadVolume(tiff_photonsFile);
  spdc_photons0 = single(spdc_photons0);
  % spdc_photons0 = spdc_photons0 / 255;
end
