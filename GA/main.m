%% ��ʼ
clear; clc;
false = 0;
true = 1;
%% ��ʼֵ + ����
N = 0; %��ʼ����������0��20��
N_total = 40; %���������
Chrom = zeros (N_total, 18);
y_max = 30;
a_max = 1.75;
a_min = 2 - a_max;
% ����Ⱥ�е�n��dʵ����(����ǵ�һ�ν�����Ⱥ��n��d��ʵ��ֵ)
Chrom_dig_n = zeros(N_total, 1);
Chrom_dig_d = zeros(N_total, 1);
Chrom_dig_h = zeros(N_total, 1);
data_need = zeros(N_total, 4);
data_loc = zeros(y_max,4);
while (N < N_total)
    var = crtbp(1, 14);
    var_str = num2str(var); %�Ѷ����Ʊ��ַ���
    var_str = strrep(var_str, ' ', ''); %ȥ���ַ����ڿո�
    % var_str = deblank(var_str);        %ȥ����β�Ķ���ո�
    var_str_n = var_str(1:6);
    var_str_d = var_str(7:14);
    var_str_h = var_str(15:18);%ȡ��n,d��h��ֵ
    % var_dig_n = str2num(var_str_n); var_dig_d = str2num(var_str_n); %���ַ����������
    var_10_n = bin2dec (var_str_n);
    var_10_d = bin2dec (var_str_d);
    var_10_h = bin2dec (var_str_h);

    if (var_10_n * var_10_d) < 390 && (var_10_n * var_10_d) > 0 && (var_10_n <= 50) && (var_10_n >= 3) &&  var_10_h > 10 && var_10_h  < 22
        N = N + 1;
        Chrom_dig_n(N, 1) = var_10_n;
        Chrom_dig_d(N, 1) = var_10_d;
        Chrom_dig_h(N, 1) = var_10_h;

        for i = 1:18
            Chrom(N, i) = var(1, i);
        end

    end

end

%% ��ʼ���㣨������
for e = 1:N_total

    %% ΢���߶�
    a = [];

        d = Chrom_dig_d(e, 1);
        n = Chrom_dig_n(e, 1);
        h = Chrom_dig_h(e, 1);
        ys = (390 - (n - 1) * d) / n;
    n1 = 26;
    n2 = 3;
    n3 = 3;

    zstandard = 21;

    PointLocation = 1;

    aaa = 1;

    for points = 1:4 * n - 1

        if points == 1
            x = 540;
            y = ys / 2;
            z = 0;
            condition = 1;
        end

        if points == 2 % x±�?0
            x = 0;
            y = y;
            z = 0;
            condition = 2;
        end

        if mod(points + 2, 4) == 1 && points > 2
            x = 0;
            y = y;
            z = zstandard;
            condition = 1;
        elseif mod(points + 2, 4) == 2 && points > 2
            x = 540;
            y = y;
            z = zstandard;

            if mod(n, 2) == 1

                if aaa <= (n - 1) / 2
                    %ddd = d / 2 + (2 * aaa - 2) / (n - 3) * d;
                    ddd = (a_min*d)+((2*aaa-2)*(a_max-a_min)*d)/(n-3);
                else
                    %ddd = d / 2 + 2 * (n - aaa - 1) / (n - 3) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-3);
                end

            else

                if aaa <= n / 2
                    %ddd = d / 2 + (aaa -1) * 2 / (n - 2) * d;
                    ddd = (a_min*d) + ((2*aaa-2)*(a_max-a_min)*d)/(n-2);
                else
                    %ddd = d / 2 + (n - aaa - 1) * 2 / (n - 2) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-2);
                end

            end

            aaa = aaa + 1;
            condition = 3;
        elseif mod(points + 2, 4) == 3 && points > 2
            x = 540;

            y = y + ddd + ys;
            z = 0;
            condition = 1;
        elseif mod(points + 2, 4) == 0 && points > 2
            x = 0;
            y = y;
            z = 0;
            condition = 2;
        end

        if condition == 1 % x��
            initalp = 0;

            for nn = PointLocation:PointLocation + n1 + 1

                if nn == PointLocation + n1 + 1

                    if x == 0
                        a(1, nn) = 540;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    else
                        a(1, nn) = 0;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    end

                    break
                end

                if x == 0
                    a(1, nn) = 540 / (n1 + 1) * initalp;
                    initalp = initalp + 1;
                    a(2, nn) = y;
                    a(3, nn) = z;
                else
                    a(1, nn) = 540 - 540 / (n1 + 1) * initalp;
                    initalp = initalp + 1;
                    a(2, nn) = y;
                    a(3, nn) = z;
                end

            end

            PointLocation = PointLocation + n1 + 1;
        end

        if condition == 2 % z��
            initalp = 0;

            for nn = PointLocation:PointLocation + n2 + 1

                if nn == PointLocation + n2 + 1

                    if z == 0
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = zstandard;
                    else
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = 0;
                    end

                    break
                end

                if z == 0
                    a(1, nn) = x;
                    a(2, nn) = y;
                    a(3, nn) = zstandard / (n2 + 1) * initalp;
                    initalp = initalp + 1;
                else
                    a(1, nn) = x;
                    a(2, nn) = y;
                    a(3, nn) = zstandard - zstandard / (n2 + 1) * initalp;
                    initalp = initalp + 1;
                end

            end

            PointLocation = PointLocation + n2 + 1;
        end

        if condition == 3
            initalp = 0;

            for nn = PointLocation:PointLocation + n3 + 1

                if nn == PointLocation + n3 + 1
                    a(1, nn) = x;
                    a(2, nn) = y + ddd + ys;
                    a(3, nn) = 0;
                    break
                end

                a(1, nn) = x;
                a(2, nn) = y + (ddd + ys) / (n3 + 1) * initalp;
                a(3, nn) = zstandard - zstandard / (n3 + 1) * initalp;
                initalp = initalp + 1;

            end

            PointLocation = PointLocation + n3 + 1;
        end

    end

    z_high = 316;
    aaa = 1;

    for points = 1:4 * n - 1

        if points == 1
            x = 540;
            y = y;
            z = z_high;
            condition = 1;
        end

        if points == 2 % x±�?0
            x = 0;
            y = y;
            z = z_high;
            condition = 2;
        end

        if mod(points + 2, 4) == 1 && points > 2
            x = 0;
            y = y;
            z = z_high + zstandard;
            condition = 1;
        elseif mod(points + 2, 4) == 2 && points > 2
            x = 540;
            y = y;
            z = z_high + zstandard;

            if mod(n, 2) == 1

                if aaa <= (n - 1) / 2
                    %ddd = d / 2 + (2 * aaa - 2) / (n - 3) * d;
                    ddd = (a_min*d)+((2*aaa-2)*(a_max-a_min)*d)/(n-3);
                else
                    %ddd = d / 2 + 2 * (n - aaa - 1) / (n - 3) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-3);
                end

            else

                if aaa <= n / 2
                    %ddd = d / 2 + (aaa -1) * 2 / (n - 2) * d;
                    ddd = (a_min*d) + ((2*aaa-2)*(a_max-a_min)*d)/(n-2);
                else
                    %ddd = d / 2 + (n - aaa - 1) * 2 / (n - 2) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-2);
                end

            end

            aaa = aaa + 1;
            condition = 3;
        elseif mod(points + 2, 4) == 3 && points > 2
            x = 540;

            y = y - ddd - ys;
            z = z_high;
            condition = 1;
        elseif mod(points + 2, 4) == 0 && points > 2
            x = 0;
            y = y;
            z = z_high;
            condition = 2;
        end

        if condition == 1 % x��
            initalp = 0;

            for nn = PointLocation:PointLocation + n1 + 1

                if nn == PointLocation + n1 + 1

                    if x == 0
                        a(1, nn) = 540;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    else
                        a(1, nn) = 0;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    end

                    break
                end

                if x == 0
                    a(1, nn) = 540 / (n1 + 1) * initalp;
                    initalp = initalp + 1;
                    a(2, nn) = y;
                    a(3, nn) = z;
                else
                    a(1, nn) = 540 - 540 / (n1 + 1) * initalp;
                    initalp = initalp + 1;
                    a(2, nn) = y;
                    a(3, nn) = z;
                end

            end

            PointLocation = PointLocation + n1 + 1;
        end

        if condition == 2 % z��
            initalp = 0;

            for nn = PointLocation:PointLocation + n2 + 1

                if nn == PointLocation + n2 + 1

                    if z == z_high
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = z_high + zstandard;
                    else
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = z_high;
                    end

                    break
                end

                if z == z_high
                    a(1, nn) = x;
                    a(2, nn) = y;
                    a(3, nn) = z_high + zstandard / (n2 + 1) * initalp;
                    initalp = initalp + 1;
                else
                    a(1, nn) = x;
                    a(2, nn) = y;
                    a(3, nn) = z_high + zstandard - zstandard / (n2 + 1) * initalp;
                    initalp = initalp + 1;
                end

            end

            PointLocation = PointLocation + n2 + 1;
        end

        if condition == 3
            initalp = 0;

            for nn = PointLocation:PointLocation + n3 + 1

                if nn == PointLocation + n3 + 1
                    a(1, nn) = x;
                    a(2, nn) = y - ddd - ys;
                    a(3, nn) = z_high;
                    break
                end

                a(1, nn) = x;
                a(2, nn) = y - (ddd + ys) / (n3 + 1) * initalp;
                a(3, nn) = z_high + zstandard - zstandard / (n3 + 1) * initalp;
                initalp = initalp + 1;

            end

            PointLocation = PointLocation + n3 + 1;
        end

    end

    %% ȡ�۲��
    %O.= [270;195;169]
    mun = 180 * 19 + 2;
    star_1 = 0;
    star_2 = 0;
    z_o = 169;
    location = zeros(3, mun);

    for h = -90:10:90
        location(3, star_1 + 1:star_1 + 180) = z_o + h;
        r = ((100)^2 - (abs(h))^2)^(1/2);

        for i = 0:1:179
            star_2 = star_2 + 1;
            location(1, star_2) = 270 + r * cosd(i * 2);
            location(2, star_2) = 195 + r * sind(i * 2);
        end

        star_1 = star_1 + 180;
    end

    location(1, mun - 1) = 270;
    location(2, mun - 1) = 195;
    location(3, mun - 1) = 169;
    location(1, mun) = 270;
    location(2, mun) = 195;
    location(3, mun) = 69;
    %% ȡ�����ȶ�
    MagData = MagneticField(a, location);
    data_need(e, 1) = max(max(MagData(2, :)));
    data_need(e, 2) = min(min(MagData(2, :)));
    data_need(e, 3) = (data_need(e, 1) - data_need(e, 2)) / ((data_need(e, 1) + data_need(e, 2)) / 2);
    fprintf('Initial Child %d is done  ', e);
    T = clock;
    fprintf('Time: %d:%d:%.0f \n', T(4), T(5), T(6));
end

%% ��ʼ�Ӵ�����
fprintf('Time: %d:%d:%.0f : Begin to calculate\n', T(4), T(5), T(6));

for yr = 1:y_max
    fprintf('Time: %d:%d:%.0f : Begin to calculate Year yr:%d\n', T(4), T(5), T(6), yr);
    %% ���븸���������Ӵ���Ⱥ
    %��ز���
    GGAP = 0.9; %����
    year = num2str(yr);
    Pm_1 = 0.1; %�������
    Pm_2 = 0.1;
    Pm_3 = 0.1;
    k_1 = 3; %n����Ľ����
    k_2 = 4; %d����Ľ����
    k_3 = 3; %d����Ľ����
    n_max = 100;
    d_max = 200;
    SUBPOP = 1; %��P88��
    InsOpt = [0, 1]; %�Ӵ������ʽ ��P88��
    %for i = 1:N_total
    %    data_need(i, 3) = (data_need(i, 1) - data_need(i, 2)) / ((data_need(i, 1) + data_need(i, 2)) / 2);
    %end
    [min_1, loc] = min(data_need(:, 3));
    FintV = ranking(data_need(:, 3)); %������Ӧ��ֵ
    SelCh = select('sus', Chrom, FintV); %ѡ��
    %[m,n]=size(SelCh);
    SelCh_str = num2str(SelCh);
    SelCh_str = char(strrep(cellstr(SelCh_str), ' ', '')); %ȥ���ַ����ڿո�
    for i = 1:N_total
        SelCh_str_n(i, 1:6) = SelCh_str(i, 1:6);
        SelCh_str_d(i, 1:8) = SelCh_str(i, 7:14);
        SelCh_str_h(i, 1:4) = SelCh_str(i, 15:18);
    end
    %ȥ�ո�
    %SelCh_dig_n = zeros (x,1); SelCh_dig_d = zeros (x,1);
    SelCh_dig_n = zeros(N_total, 6); SelCh_dig_d = zeros(N_total, 8);SelCh_dig_h = zeros(N_total, 4);

    for i = 1:N_total

        for n_1 = 1:6
            SelCh_dig_n(i, n_1) = str2double(SelCh_str_n(i, n_1));
            %SelCh_str_n_new(i,:) = dec2bin(SelCh_dig_n(i,:));
        end

        for m_1 = 1:8
            SelCh_dig_d(i, m_1) = str2double(SelCh_str_d(i, m_1));
            %SelCh_str_d_new(i,:) = dec2bin(SelCh_dig_d(i,:));
        end
        
        for h_1 = 1:4
            SelCh_dig_h(i, h_1) = str2double(SelCh_str_h(i, h_1));
            %SelCh_str_d_new(i,:) = dec2bin(SelCh_dig_d(i,:));
        end

    end

    %SelCh_str_d = char(strrep(cellstr(SelCh_str_d),' ',''));  %ȥ���ַ����ڿո�
    %recombin part
    %SelCh_n = recombin('xovsp',SelCh_str_n);
    %SelCh_d = recombin('xovsp',SelCh_str_d);
    %�������
    SelCh_n = xovsp (SelCh_dig_n, k_1);
    SelCh_d = xovsp (SelCh_dig_d, k_2);
    SelCh_h = xovsp (SelCh_dig_h, k_3);
    %����
    SelCh_n = mut(SelCh_n, Pm_1);
    SelCh_d = mut(SelCh_d, Pm_2);
    SelCh_h = mut(SelCh_h, Pm_3);
    SelCh_new = [SelCh_n, SelCh_d,SelCh_h]; %�ѻ���Ƭ���������
    %����
    %SelCh_new = mutate('mut',SelCh_new);
    %Chrom = reins(Chrom, SelCh, SUBPOP, InsOpt);
    Chrom = SelCh_new;
    %��������
    data_need = zeros(N_total, 3);
    Chrom_dig_n = zeros(N_total, 1);
    Chrom_dig_d = zeros(N_total, 1);
    Chrom_dig_h = zeros(N_total, 1);
    Chrom_str = num2str(Chrom);
    Chrom_str = char(strrep(cellstr(Chrom_str), ' ', '')); %ȥ���ַ����ڿո�

    for i = 1:N_total
        Chrom_10_n = Chrom_str(i, 1:6);
        Chrom_dig_n(i, 1) = bin2dec (Chrom_10_n);
        Chrom_10_d = Chrom_str(i, 7:14);
        Chrom_dig_d(i, 1) = bin2dec (Chrom_10_d);
        Chrom_10_h = Chrom_str(i,15:18);
        Chrom_dig_h(i, 1) = bin2dec (Chrom_10_h);
    end

    for i = 1:N_total

        if Chrom_dig_n(i, 1) <= 3 || Chrom_dig_n(i, 1) > n_max || (Chrom_dig_n(i, 1) - 1) * Chrom_dig_d(i, 1) > 390 || Chrom_dig_d(i, 1) <= 0 || Chrom_dig_d(i, 1) > d_max || Chrom_dig_h(i, 1) < 10|| Chrom_dig_h(i, 1) > 22 
            N_1 = 0;

            while N_1 <= 0
                var_new = crtbp(1, 14);
                var_new_str = num2str(var_new); %�Ѷ����Ʊ��ַ���
                var_new_str = strrep(var_new_str, ' ', ''); %ȥ���ַ����ڿո�
                var_new_str_n = var_new_str(1:6);
                var_new_str_d = var_new_str(7:14); 
                var_new_str_h = var_new_str(15:18); %ȡ��n��d��ֵ
                var_new_10_n = bin2dec (var_new_str_n);
                var_new_10_d = bin2dec (var_new_str_d);
                var_new_10_h = bin2dec (var_new_str_h);

                if (var_new_10_n * var_new_10_d) < 390 && (var_new_10_n * var_new_10_d) > 0 && (var_new_10_n <= 50) && (var_new_10_n >= 3)
                    N_1 = N_1 + 1;
                    Chrom_dig_n(i, 1) = var_new_10_n;
                    Chrom_dig_d(i, 1) = var_new_10_d;
                    Chrom_dig_h(i, 1) = var_new_10_h;

                    for j = 1:18
                        Chrom(i, j) = var_new(1, j);
                    end

                end

            end

            %            Chrom_dig_n(i, 1) = bin2dec (SelCh_str_n(loc,:));
            %            Chrom(i,1:6) = SelCh(loc,1:6);
            %        end
            %        if Chrom_dig_d(i, 1)<= 0 || Chrom_dig_d(i, 1)> d_max || Chrom_dig_n(i, 1)*Chrom_dig_d(i, 1) > 390
            %            Chrom_dig_d(i, 1) = bin2dec (SelCh_str_d(loc,:));
            %            Chrom(i,7:14) = SelCh(loc,7:14);
        end

    end

    %% ��ʼ�����Ӵ�
    for e = 1:N_total
        %% ΢���߶�
        %% ΢���߶�
        a = [];
        d = Chrom_dig_d(e, 1);
        n = Chrom_dig_n(e, 1);
        ys = (390 - (n - 1) * d) / n;
        n1 = 26;
        n2 = 3;
        n3 = 3;

        zstandard = 21;
        aaa = 1;
        PointLocation = 1;

        for points = 1:4 * n - 1

            if points == 1
                x = 540;
                y = ys / 2;
                z = 0;
                condition = 1;
            end

            if points == 2 % x±�?0
                x = 0;
                y = y;
                z = 0;
                condition = 2;
            end

            if mod(points + 2, 4) == 1 && points > 2
                x = 0;
                y = y;
                z = zstandard;
                condition = 1;
            elseif mod(points + 2, 4) == 2 && points > 2
                x = 540;
                y = y;
                z = zstandard;

                if mod(n, 2) == 1

                    if aaa <= (n - 1) / 2
                        %ddd = d / 2 + (2 * aaa - 2) / (n - 3) * d;
                        ddd = (a_min*d)+((2*aaa-2)*(a_max-a_min)*d)/(n-3);
                    else
                        %ddd = d / 2 + 2 * (n - aaa - 1) / (n - 3) * d;
                        ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-3);
                    end

                else

                    if aaa <= n / 2
                        %ddd = d / 2 + (aaa -1) * 2 / (n - 2) * d;
                        ddd = (a_min*d)+((2*aaa-2)*(a_max-a_min)*d)/(n-2);
                    else
                        %ddd = d / 2 + (n - aaa - 1) * 2 / (n - 2) * d;
                        ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-2);
                    end

                end

                aaa = aaa + 1;
                condition = 3;
            elseif mod(points + 2, 4) == 3 && points > 2
                x = 540;
                y = y + ddd + ys;
                z = 0;
                condition = 1;
            elseif mod(points + 2, 4) == 0 && points > 2
                x = 0;
                y = y;
                z = 0;
                condition = 2;
            end

            if condition == 1 % x��
                initalp = 0;

                for nn = PointLocation:PointLocation + n1 + 1

                    if nn == PointLocation + n1 + 1

                        if x == 0
                            a(1, nn) = 540;
                            a(2, nn) = y;
                            a(3, nn) = z;
                        else
                            a(1, nn) = 0;
                            a(2, nn) = y;
                            a(3, nn) = z;
                        end

                        break
                    end

                    if x == 0
                        a(1, nn) = 540 / (n1 + 1) * initalp;
                        initalp = initalp + 1;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    else
                        a(1, nn) = 540 - 540 / (n1 + 1) * initalp;
                        initalp = initalp + 1;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    end

                end

                PointLocation = PointLocation + n1 + 1;
            end

            if condition == 2 % z��
                initalp = 0;

                for nn = PointLocation:PointLocation + n2 + 1

                    if nn == PointLocation + n2 + 1

                        if z == 0
                            a(1, nn) = x;
                            a(2, nn) = y;
                            a(3, nn) = zstandard;
                        else
                            a(1, nn) = x;
                            a(2, nn) = y;
                            a(3, nn) = 0;
                        end

                        break
                    end

                    if z == 0
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = zstandard / (n2 + 1) * initalp;
                        initalp = initalp + 1;
                    else
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = zstandard - zstandard / (n2 + 1) * initalp;
                        initalp = initalp + 1;
                    end

                end

                PointLocation = PointLocation + n2 + 1;
            end

            if condition == 3
                initalp = 0;

                for nn = PointLocation:PointLocation + n3 + 1

                    if nn == PointLocation + n3 + 1
                        a(1, nn) = x;
                        a(2, nn) = y + ddd + ys;
                        a(3, nn) = 0;
                        break
                    end

                    a(1, nn) = x;
                    a(2, nn) = y + (ddd + ys) / (n3 + 1) * initalp;
                    a(3, nn) = zstandard - zstandard / (n3 + 1) * initalp;
                    initalp = initalp + 1;

                end

                PointLocation = PointLocation + n3 + 1;
            end

        end

        aaa = 1;
        z_high = 316;

        for points = 1:4 * n - 1

            if points == 1
                x = 540;
                y = y;
                z = z_high;
                condition = 1;
            end

            if points == 2 % x±�?0
                x = 0;
                y = y;
                z = z_high;
                condition = 2;
            end

            if mod(points + 2, 4) == 1 && points > 2
                x = 0;
                y = y;
                z = z_high + zstandard;
                condition = 1;
            elseif mod(points + 2, 4) == 2 && points > 2
                x = 540;
                y = y;
                z = z_high + zstandard;

            if mod(n, 2) == 1

                if aaa <= (n - 1) / 2
                    %ddd = d / 2 + (2 * aaa - 2) / (n - 3) * d;
                    ddd = (a_min*d)+((2*aaa-2)*(a_max-a_min)*d)/(n-3);
                else
                    %ddd = d / 2 + 2 * (n - aaa - 1) / (n - 3) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-3);
                end

            else

                if aaa <= n / 2
                    %ddd = d / 2 + (aaa -1) * 2 / (n - 2) * d;
                    ddd = (a_min*d) + ((2*aaa-2)*(a_max-a_min)*d)/(n-2);
                else
                    %ddd = d / 2 + (n - aaa - 1) * 2 / (n - 2) * d;
                    ddd = (a_min*d)+(2*(n-aaa-1)*(a_max-a_min)*d)/(n-2);
                end

            end

                aaa = aaa + 1;
                condition = 3;
            elseif mod(points + 2, 4) == 3 && points > 2
                x = 540;
                y = y - ddd - ys;
                z = z_high;
                condition = 1;
            elseif mod(points + 2, 4) == 0 && points > 2
                x = 0;
                y = y;
                z = z_high;
                condition = 2;
            end

            if condition == 1 % x��
                initalp = 0;

                for nn = PointLocation:PointLocation + n1 + 1

                    if nn == PointLocation + n1 + 1

                        if x == 0
                            a(1, nn) = 540;
                            a(2, nn) = y;
                            a(3, nn) = z;
                        else
                            a(1, nn) = 0;
                            a(2, nn) = y;
                            a(3, nn) = z;
                        end

                        break
                    end

                    if x == 0
                        a(1, nn) = 540 / (n1 + 1) * initalp;
                        initalp = initalp + 1;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    else
                        a(1, nn) = 540 - 540 / (n1 + 1) * initalp;
                        initalp = initalp + 1;
                        a(2, nn) = y;
                        a(3, nn) = z;
                    end

                end

                PointLocation = PointLocation + n1 + 1;
            end

            if condition == 2 % z��
                initalp = 0;

                for nn = PointLocation:PointLocation + n2 + 1

                    if nn == PointLocation + n2 + 1

                        if z == z_high
                            a(1, nn) = x;
                            a(2, nn) = y;
                            a(3, nn) = z_high + zstandard;
                        else
                            a(1, nn) = x;
                            a(2, nn) = y;
                            a(3, nn) = z_high;
                        end

                        break
                    end

                    if z == z_high
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = z_high + zstandard / (n2 + 1) * initalp;
                        initalp = initalp + 1;
                    else
                        a(1, nn) = x;
                        a(2, nn) = y;
                        a(3, nn) = z_high + zstandard - zstandard / (n2 + 1) * initalp;
                        initalp = initalp + 1;
                    end

                end

                PointLocation = PointLocation + n2 + 1;
            end

            if condition == 3
                initalp = 0;

                for nn = PointLocation:PointLocation + n3 + 1

                    if nn == PointLocation + n3 + 1
                        a(1, nn) = x;
                        a(2, nn) = y - ddd - ys;
                        a(3, nn) = z_high;
                        break
                    end

                    a(1, nn) = x;
                    a(2, nn) = y - (ddd + ys) / (n3 + 1) * initalp;
                    a(3, nn) = z_high + zstandard - zstandard / (n3 + 1) * initalp;
                    initalp = initalp + 1;

                end

                PointLocation = PointLocation + n3 + 1;
            end

        end

        %% ȡ�۲��
        %O.= [270;195;169]
        mun = 180 * 19 + 2;
        star_1 = 0;
        star_2 = 0;
        z_o = 169;
        location = zeros(3, mun);

        for h = -90:10:90
            location(3, star_1 + 1:star_1 + 180) = z_o + h;
            r = ((100)^2 - (abs(h))^2)^(1/2);

            for i = 0:1:179
                star_2 = star_2 + 1;
                location(1, star_2) = 270 + r * cosd(i * 2);
                location(2, star_2) = 195 + r * sind(i * 2);
            end

            star_1 = star_1 + 180;
        end

        location(1, mun - 1) = 270;
        location(2, mun - 1) = 195;
        location(3, mun - 1) = 169;
        location(1, mun) = 270;
        location(2, mun) = 195;
        location(3, mun) = 69;
        %% ȡ�����ȶ�
        MagData = MagneticField(a, location);
        data_need(e, 1) = max(max(MagData(2, :)));
        data_need(e, 2) = min(min(MagData(2, :)));
        data_need(e, 3) = (data_need(e, 1) - data_need(e, 2)) / ((data_need(e, 1) + data_need(e, 2)) / 2);
        FileName = 'Report_Year_';
        FileName = strcat(FileName, year);
        FileName = strcat(FileName, '.csv');
        T = table(Chrom_dig_n, Chrom_dig_d, data_need(:, 1), data_need(:, 2), data_need(:, 3));
        writetable(T, FileName)

        fprintf('Child %d is done!  ', e);
        T = clock;
        fprintf('Time: %d:%d:%.0f \n', T(4), T(5), T(6));
    end
        [min_2,loc_1] = min(data_need(:,3));
        data_loc(yr,1) = yr; 
        data_loc(yr,2) = Chrom_dig_n(loc_1,1); 
        data_loc(yr,3) = Chrom_dig_d(loc_1,1); 
        data_loc(yr,4) = min_2; 
        fprintf('\nYear:%d n = %d d = %d Min = %.6f\n', data_loc(yr,1),data_loc(yr,2),data_loc(yr,3),data_loc(yr,4));
        TT = table(data_loc(:,1),data_loc(:,2),data_loc(:,3),data_loc(:,4));
        writetable(TT, "��ȱ�.csv")
    fprintf('%s is done!\n', FileName);
    T = clock;
    fprintf('Time: %d:%d:%.0f \n\n\n', T(4), T(5), T(6));
end
