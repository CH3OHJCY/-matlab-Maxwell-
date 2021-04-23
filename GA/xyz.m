clc;clear;
%O.= [270;195;169]
mun = 180 * 19 + 2;
star_1 = 0;
star_2 = 0;
z_o = 169;
points = zeros(3,mun);
for d = -90:10:90
    points(3,star_1+1:star_1+180) = z_o + d;
    r = ((100)^2 - (abs(d))^2)^(1/2);
    for a = 0 : 1 : 179
        star_2 = star_2+1;
        points(1,star_2) = 270+r*cosd(a*2);
        points(2,star_2) = 195+r*sind(a*2);
    end
    star_1 = star_1 + 180;
end
points(1,mun-1) = 270;points(2,mun-1) = 195;points(3,mun-1) = 169;
points(1,mun) = 270;  points(2,mun) = 195;  points(3,mun) = 69;