%% 开始
clear; clc;
false = 0;
true = 1;
%% 创建初始种群
N = 0; %初始个体数（从0到20）
N_total = 20; %总体个体数
Chrom = zeros (N_total, 14);
% 把种群中的n和d实数化(这个是第一次交叉种群的n和d的实数值)
Chrom_dig_n = zeros(N_total, 1);
Chrom_dig_d = zeros(N_total, 1);
data_need = zeros(N_total, 3);

while (N < N_total)
    var = crtbp(1, 14);
    var_str = num2str(var); %把二进制变字符串
    var_str = strrep(var_str, ' ', ''); %去除字符串内空格
    % var_str = deblank(var_str);        %去除首尾的多余空格
    var_str_n = var_str(1:6);
    var_str_d = var_str(7:14); %取出n和d的值
    % var_dig_n = str2num(var_str_n); var_dig_d = str2num(var_str_n); %把字符串变二进制
    var_10_n = bin2dec (var_str_n);
    var_10_d = bin2dec (var_str_d);

    if (var_10_n * var_10_d) < 390 && (var_10_n * var_10_d) > 0 && (var_10_n <= 10) &&  (var_10_n >= 3)
        N = N + 1;
        Chrom_dig_n(N, 1) = var_10_n;
        Chrom_dig_d(N, 1) = var_10_d;

        for i = 1:14
            Chrom(N, i) = var(1, i);
        end

    end

end

% 输出.mat文件
save data_fc  Chrom_dig_n Chrom_dig_d data_need


% 输出.mat文件
save data_fc  Chrom_dig_n Chrom_dig_d data_need Chrom N_total;