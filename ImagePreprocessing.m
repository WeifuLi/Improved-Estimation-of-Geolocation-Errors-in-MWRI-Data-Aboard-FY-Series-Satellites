function [T5,T6]=ImagePreprocessing(LandSeaMask)
T1=LandSeaMask==1; 
T2=LandSeaMask==2;
T1=T1|T2; %将陆地填充淡水区域
SE=strel('disk',1);
T11=imdilate(T1,SE);
T12=bwareaopen(T11,4000);
T3=LandSeaMask==3;
T31=imdilate(T3,SE);
T32=bwareaopen(T31,4000);
T=T12&T32; %figure;imshow(T,[]);
T5=LandSeaMask==5;
T6=T;