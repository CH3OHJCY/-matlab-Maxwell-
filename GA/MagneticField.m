function [B] = MagneticField(cpath,tpos,I)
%MAGNETICFIELD ����ͨ�絼����ָ�������ڲ����Ĵų��ֲ�
%   [B] = MagneticField(cpath,tpos);
%   [B] = MagneticField(cpath,tpos,I);
%
%   ���������
%           cpath	����·���������ݣ�1,2,3�зֱ��Ӧx,y,z����
%           tpos	����Ŀ����������ݣ�1,2,3�зֱ��Ӧx,y,z����
%           I       ������ͨ�������С
%   ���������
%           B   �ų�ֵ
%
% Author: He Yucheng
% date: 2020-07-07
% version: v1.0
% Email: heyucheng@cqu.edu.cn

if nargin==2        % �ж������������
    I = 1;
elseif nargin<=1
    error('Please input correct parameters!\nUsage: MagnticField(cpath,tpos,I)',nargin);
end

u0 = pi*4e-7;		% ��մŵ���
npos = size(tpos);      % 
ncp = length(cpath);        % ����·������
B = zeros(3,npos(2));       % ��ʼ��B

for m=1:npos(2)         
    for k=1:(ncp-1)         
        dl(:,k) = cpath(:,k+1)-cpath(:,k);    % dl�߶λ���΢Ԫ
        dlpos(:,k) = (cpath(:,k+1)+cpath(:,k))/2;   % dl�߶���������
        r = dlpos(:,k) - tpos(:,m);
        dB = u0.*I./(4*pi)*cross(dl(:,k),r)./magnitude(r)^3;    % Biot-Savart law
        B(:,m) = B(:,m)+dB;
    end
end
end

