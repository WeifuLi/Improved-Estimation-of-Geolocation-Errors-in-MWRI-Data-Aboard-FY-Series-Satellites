% ȥë��
function [BwRegionValue]=BwRemoveSpur(BwRegionValue)
% �ҷֲ��
Bwbranchpoints=bwmorph(BwRegionValue,'branchpoints');
if sum(Bwbranchpoints(:))==1
    BwRegionValue1=BwRegionValue;
    % �Էֲ��Ϊ�ڵ�Ͽ����������ɸ���ͨ��
    Bwbranchpoints1=bwmorph(Bwbranchpoints,'dilate');
    BwRegionValue1(Bwbranchpoints1)=0;
    [Bw,~]=bwlabel(BwRegionValue1,8);
    State=regionprops(Bw,'Area');
    Areas=[State.Area];
    [~,IndM]=min(Areas);
    BwRegionValue(Bw==IndM)=0;
end

if sum(Bwbranchpoints(:))>1
    BwRegionValue1=BwRegionValue;
    % �Էֲ��Ϊ�ڵ�Ͽ����������ɸ���ͨ��
    Bwbranchpoints1=bwmorph(Bwbranchpoints,'dilate');
    BwRegionValue1(Bwbranchpoints1)=0;
    [Bw,Num]=bwlabel(BwRegionValue1,8);
    State=regionprops(Bw,'Area');
    Areas=[State.Area];
    Reserved=[];
    for Ii=1:Num
        TemNw=Bw==Ii;
        TemNw=bwmorph(TemNw,'dilate',2);
        TemNw=Bwbranchpoints&TemNw;
        % ���������зֲ�����ͨ��
        if sum(TemNw(:))>1
            Reserved=[Reserved Ii];
        end
    end
    [IxB,IyB]=find(Bwbranchpoints==1);
    for Ii=1:length(IxB)
        TemIxB=IxB(Ii);TemIyB=IyB(Ii);
        if TemIxB<4||TemIxB>size(Bwbranchpoints,1)-3||TemIyB<4||TemIyB>size(Bwbranchpoints,2)-3
            continue;
        end
        TemR=BwRegionValue(TemIxB-3:TemIxB+3,TemIyB-3:TemIyB+3);
        TemR(3:5,3:5)=0;
        BwV=Bw(TemIxB-3:TemIxB+3,TemIyB-3:TemIyB+3);
        V=BwV(TemR==1);
        V1=unique(V);
        V1=setdiff(V1,[0 Reserved]);
        % �������˵㴦�����ͨ��  ��Ҫ�㹻��
        if length(V1)>1
            [MaxV1,IndM]=max(Areas(V1));
            if MaxV1>30
                Reserved=[Reserved V1(IndM)];
            end
        end
    end
    Removed=setdiff(1:Num,Reserved);
    for Ii=1:length(Removed)
        BwInd=Bw==Removed(Ii);
        se=strel('diamond',1);
        BwInd=imdilate(BwInd,se);
        BwRegionValue(BwInd)=0;
    end
end
end