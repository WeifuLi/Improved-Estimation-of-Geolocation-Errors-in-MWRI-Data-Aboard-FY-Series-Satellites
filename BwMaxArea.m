function [Image]=BwMaxArea(Image)
% �ҵ�����һ����ͨ��
[Bw,Num]=bwlabel(Image);
if Num>1
    State=regionprops(Bw,'Area');
    Areas=[State.Area];
    [~,IndM]=max(Areas);
    Image=Bw==IndM;
end
end