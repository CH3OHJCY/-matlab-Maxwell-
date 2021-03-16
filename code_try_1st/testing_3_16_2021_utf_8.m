%%遗传算法第一次尝试
% 2021/3/16 贾淳宇
clc;clear;clc;
%% 初始种群建立
N = 0;  %初始个体数（从0到20）
Chrom = zeros (20,14);
while (N<20)
    var = crtbp(1,14);
    var_str = num2str(var);  %把二进制变字符串 
    var_str=strrep(var_str,' ','');      %去除字符串内空格
% var_str = deblank(var_str);        %去除首尾的多余空格
    var_str_n = var_str(1:6); var_str_d = var_str(7:14); %取出n和d的值
% var_dig_n = str2num(var_str_n); var_dig_d = str2num(var_str_n); %把字符串变二进制
    var_10_n = bin2dec (var_str_n); var_10_d = bin2dec (var_str_d);
    if (var_10_n * var_10_d) < 390
        N = N+1;
        for i = 1:14
            Chrom(N,i) = var(1,i);
        end
    end
end

%% 输入 Maxwell 计算


%% Maxwell csv文件读取
%data = 'C:\Users\ch3oh\Desktop\marker_3_16_2.csv';
data = csvread('marker_3_16_2.csv',1,1); 
%1.M = csvread('filename') 
%2.M = csvread('filename', row, col) 行号、列号以0开始计数。也就是说，row=0, col=0表示从文件中第一个数开始读
%3.M = csvread('filename', row, col, range) range = [R1 C1 R2 C2]，这里（R1，C1）是读取区域的左上角，（R2，C2）是读取区域的右下角

%% 开始迭代
%相关参数
MAXGEN = 25;  %最大迭代数
GGAP = 0.9;    %代沟 
years = 0;     %迭代次数
Pm =0.1;      %变异概率
k_1 = 3;      %n基因的交叉点
k_2 = 4;      %d基因的交叉点
SUBPOP = 1;   %（P88）
InsOpt = [0,1];   %子代替代方式 （P88）
%计算不不均匀度
x = 20; %行数
y = 1; %交叉概率
Nonun = zeros (x,1);
for i = 1:x
    Nonun(i,1) = max(data(i,:)) - min (data(i,:)); 
end
FintV = ranking(-Nonun);         %分配适应度值
SelCh = select('sus',Chrom,FintV);   %选择
[m,n]=size(SelCh);
SelCh_str = num2str(SelCh);
SelCh_str = char(strrep(cellstr(SelCh_str),' ',''));  %去除字符串内空格
for i = 1:x
    SelCh_str_n(i,1:6) = SelCh_str(i,1:6);
    SelCh_str_d(i,1:8) = SelCh_str(i,7:14);
end

%去空格
%SelCh_dig_n = zeros (x,1); SelCh_dig_d = zeros (x,1);

SelCh_dig_n = zeros(20,6);  SelCh_dig_d = zeros(20,8);
for i = 1:x
    for n = 1:6
        SelCh_dig_n(i,n) = str2double(SelCh_str_n(i,n));
        %SelCh_str_n_new(i,:) = dec2bin(SelCh_dig_n(i,:));
    end
    for m = 1:8
        SelCh_dig_d(i,m) = str2double(SelCh_str_d(i,m));
        %SelCh_str_d_new(i,:) = dec2bin(SelCh_dig_d(i,:));
    end
end

%SelCh_str_d = char(strrep(cellstr(SelCh_str_d),' ',''));  %去除字符串内空格
%recombin part 
%SelCh_n = recombin('xovsp',SelCh_str_n);
%SelCh_d = recombin('xovsp',SelCh_str_d);

%交叉基因
SelCh_n = xovsprs (SelCh_dig_n,k_1);
SelCh_d = xovsprs (SelCh_dig_d,k_2);
SelCh_new =[SelCh_n,SelCh_d];         %把基因片段组合起来
%变异
SelCh_new = mut(SelCh_new,Pm);
Chrom = reins(Chrom,SelCh,SUBPOP,InsOpt);