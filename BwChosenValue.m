function [Mask]=BwChosenValue(Bw,BwValue)
% ѡȡBwͼ��ĳ����ֵ������
if isempty(BwValue)==1
    Mask=Bw>0;
else
    Mask=zeros(size(Bw));
    for i=1:length(BwValue)
        TemValue=BwValue(i);
        Mask(Bw==TemValue)=1;
    end
end
