function [Result]=MaxValueinVertical(BwRegionValue,D2,D)
% BwRegionValue 表示提取的线
% D2 表示通过形态学处理提取的线
% D  表示原始的值
Result=zeros(size(D));
Se=ones(21,21);
BwRegionValue1=conv2(BwRegionValue,Se,'same');
BwRegionValue=BwRegionValue&(BwRegionValue1>20);
[Indx,Indy]=find(BwRegionValue==1);
r=2;
for i=1:length(Indx)
    Ix=Indx(i);
    Iy=Indy(i);
    if Ix<r+1||Iy<r+1||Ix>size(BwRegionValue,1)-r||Iy>size(BwRegionValue,2)-r
        continue;
    end
    Ix1=Ix-r:Ix+r;
    Iy1=Iy-r:Iy+r;
    BwRegion=BwRegionValue(Ix1,Iy1);
    if sum(BwRegion(:))<5
        continue;
    end
    Mask=zeros(size(D));
    Mask(Ix,Iy)=1;
    [Ix2,Iy2]=find(BwRegion==1);
    [P]=pca([Ix2 Iy2]);
    Theta=180*atan(P(2,1)/P(1,1))/pi;
    SE = strel('line', 5, Theta);
    Mask=imdilate(Mask,SE);
    
    Mask=Mask&D2;
    IndM=find(Mask==1);
    TemBT=D(IndM);
    [~,Ind]=max(TemBT);
    
    Result(IndM(Ind))=1;
end
end