function [Image]=BwMaxArea(Image)
% 找到最大的一个连通域
[Bw,Num]=bwlabel(Image);
if Num>1
    State=regionprops(Bw,'Area');
    Areas=[State.Area];
    [~,IndM]=max(Areas);
    Image=Bw==IndM;
end
end