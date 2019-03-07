%% 根据定位的误差调整经纬度的值
clc;clear;close all;
RegionName='阿拉伯海岸';
AscendPath=['C:\Users\huj\Desktop\曲面求拐点\曲面求拐点\数据\' RegionName '\'];
OBCPath=['C:\Users\huj\Desktop\曲面求拐点\曲面求拐点\数据\' RegionName 'OBC\'];
TifPath='C:\Users\huj\Desktop\曲面求拐点\曲面求拐点\数据\';
ErrorPath=['修正前误差' RegionName  '\'];
mkdir(ErrorPath);
Dir=dir([AscendPath, '*.hdf']);
% load gshhs_land_f2;
load TruthCoastline;
TotalErrors=[];
TotalPixelErrors=[];
for i=1:length(Dir)
    data=hdf5info([AscendPath Dir(i).name]);
    LandSeaMask=hdf5read(data.GroupHierarchy.Groups(1).Datasets(4));
    Latitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(1)));
    Longitude=double(hdf5read(data.GroupHierarchy.Groups(2).Datasets(2)));
    ObserveData=hdf5read(data.GroupHierarchy.Groups(1).Datasets(2));
    ObserveData=double(ObserveData(:,:,10));
    ObserveData=double(ObserveData)*0.01+327.68;%转变为正值
    
    Filename=data.Filename;
    DATA.Latitude=Latitude;
    DATA.Longitude=Longitude;
    DATA.ObserveData=ObserveData;
    DATA.LandSeaMask=LandSeaMask;
    %     DATA.TruthLatitude=TruthLatitude;
    %     DATA.TruthLongitude=TruthLongitude;
    DATA.TruthCoastline=TruthCoastline;
    DATA.Filename=Filename;
    Method='Before';
    [HDF, Errors,PixelErrors]=ImprovedGER(DATA,RegionName,Method);
    TotalErrors=[TotalErrors;Errors];
    TotalPixelErrors=[TotalPixelErrors; PixelErrors];
    HDF.Closest_latitudes=HDF.Closest_latitudes+0.26;        %最近点经度
    HDF.Closest_longitudes=HDF.Closest_longitudes-0.02;      %最近点纬度
    hdf_create(HDF);
    Dst_Path=[AscendPath,'需要修正的文件\'];
    mkdir(Dst_Path); 
    copyfile([AscendPath Filename],Dst_Path);%复制当前文件
    movefile([Filename(1:end-4)  '1.HDF'],Dst_Path);%复制生成的文件
    copyfile([OBCPath Filename(1:end-13) '_OBCXX_MS.HDF'],Dst_Path);%复制OBC文件
    figure;scatter(Errors(:,2),Errors(:,1),'b.');
    hold on;scatter(mean(Errors(:,2)),mean(Errors(:,1)),'r*','LineWidth',1.5);
    axis([-0.16 0.16 -0.16 0.16]);
    xlabel('Cross-track');ylabel('Along-track');
    grid on;
    title(['点的个数' num2str(size(Errors,1)) '; 修正前: lat： ' num2str(mean(Errors(:,1)),2), '  Lon: ' num2str(mean(Errors(:,2)),2) ]);
    saveas(gcf,[ErrorPath, RegionName Dir(i).name(1:end-4) '修正前.tif'  ]);
    save ([ErrorPath, Dir(i).name(1:end-4) '修正前.mat'], 'Errors','PixelErrors') ;
end
save ([ErrorPath, 'TotalErrors修正前.mat'], 'TotalErrors','TotalPixelErrors') ;
figure;scatter(TotalErrors(:,2),TotalErrors(:,1),'b.');
hold on;scatter(mean(TotalErrors(:,2)),mean(TotalErrors(:,1)),'r*','LineWidth',1.5);
axis([-0.16 0.16 -0.16 0.16]);     xlabel('Cross-track');ylabel('Along-track');
grid on;

% 消除孤立点
TotalErrors=TotalErrors(:,1:2);
Mean=mean(TotalErrors);
Delta=sqrt(sum((TotalErrors-Mean).^2,2));
Ind=find(Delta<2*mean(Delta));
TotalErrorsN=TotalErrors(Ind,:);
TotalPixelErrors=TotalPixelErrors(Ind,:);
% H=dendrogram(Z2); 
figure;scatter(TotalErrorsN(:,2),TotalErrorsN(:,1),'b.');
hold on;scatter(mean(TotalErrorsN(:,2)),mean(TotalErrorsN(:,1)),'r*','LineWidth',1.5);
xlabel('Cross-track');ylabel('Along-track');
axis([-0.16 0.16 -0.16 0.16]);
grid on;
save ([RegionName 'TotalErrors.mat'], 'TotalErrors','TotalPixelErrors');

