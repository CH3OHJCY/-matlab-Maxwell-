%% 开始
clear; clc;
false = 0;
true = 1;
%% 开始迭代
%相关参数
GGAP = 0.9;    %代沟 
years = 1;     %迭代次数
Pm =0.1;      %变异概率
k_1 = 3;      %n基因的交叉点
k_2 = 4;      %d基因的交叉点
SUBPOP = 1;   %（P88）
InsOpt = [0,1];   %子代替代方式 （P88）
load data_fc;
Nonun = zeros (N_total,1);
for i = 1:N_total
    Nonun(i,1) = (data_need(i,1)-data_need(i,2))/((data_need(i,1)+data_need(i,2))/2); 
end
FintV = ranking(-Nonun);         %分配适应度值
SelCh = select('sus',Chrom,FintV);   %选择
%[m,n]=size(SelCh);
SelCh_str = num2str(SelCh);
SelCh_str = char(strrep(cellstr(SelCh_str),' ',''));  %去除字符串内空格

for i = 1:N_total
    SelCh_str_n(i,1:6) = SelCh_str(i,1:6);
    SelCh_str_d(i,1:8) = SelCh_str(i,7:14);
end

%去空格
%SelCh_dig_n = zeros (x,1); SelCh_dig_d = zeros (x,1);

SelCh_dig_n = zeros(N_total,6);  SelCh_dig_d = zeros(N_total,8);
for i = 1:N_total
    for n_1 = 1:6
        SelCh_dig_n(i,n_1) = str2double(SelCh_str_n(i,n_1));
        %SelCh_str_n_new(i,:) = dec2bin(SelCh_dig_n(i,:));
    end
    for m_1 = 1:8
        SelCh_dig_d(i,m_1) = str2double(SelCh_str_d(i,m_1));
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
%数据重置
data_need = zeros(N_total, 3);
Chrom_dig_n = zeros(N_total, 1);
Chrom_dig_d = zeros(N_total, 1);
Chrom_str = num2str(Chrom);
Chrom_str = char(strrep(cellstr(Chrom_str),' ',''));  %去除字符串内空格
for i = 1:N_total
    Chrom_10_n = Chrom_str(i,1:6);
    Chrom_dig_n(i,1) = bin2dec (Chrom_10_n);
    Chrom_10_d = Chrom_str(i,7:14);
    Chrom_dig_d(i,1) = bin2dec (Chrom_10_d);
end
save data_fc  Chrom_dig_n Chrom_dig_d data_need Chrom N_total;