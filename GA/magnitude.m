function [mag] = magnitude(vector)
%MAGNITUDE ����������ģֵ
%   [mag] = magnitude(vector)
%
%   vector     ����
%   mag        ģֵ
%
% Author: Yucheng He
% date: 2020-07-07
% version: v1.0
% Email: heyucheng@cqu.edu.cn

    square = vector.*vector;
    mag = sqrt(sum(square));
end

