clc;clear;close all;
RegionName='利比亚';
AscendPath=[RegionName '\'];
TifPath='.\';
Dir=dir([AscendPath, '*.hdf']);
load TruthCoastline1;
TruthCoastline=TruthCoastline1;
load TruthCoastline2;
TruthCoastline=[TruthCoastline; TruthCoastline2];
TotalErrors=[];
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
    DATA.TruthCoastline=TruthCoastline;
    DATA.Filename=Filename;
    Method='Before';
    [HDF, Errors]=ImprovedGER_Limited(DATA,RegionName,Method);
    figure;scatter(Errors(:,2),Errors(:,1),'b.');
    hold on;scatter(mean(Errors(:,2)),mean(Errors(:,1)),'r*','LineWidth',1.5);
    axis([-0.16 0.16 -0.16 0.16]);
    xlabel('Cross-track');ylabel('Along-track');
    grid on;
    title(['点的个数' num2str(size(Errors,1)) '; 修正前: lat： ' num2str(mean(Errors(:,1)),2), '  Lon: ' num2str(mean(Errors(:,2)),2) ]);
    saveas(gcf,[TifPath, RegionName Dir(i).name(1:end-4) '修正前.tif'  ]);
    save ([AscendPath, Dir(i).name(1:end-4) '.mat'], 'Errors') ;
    TotalErrors=[TotalErrors; Errors i*ones(size(Errors,1),1)];
end
save ([AscendPath, 'TotalErrors.mat'], 'TotalErrors') ;
figure;scatter(TotalErrors(:,2),TotalErrors(:,1),'b.');
hold on;scatter(mean(TotalErrors(:,2)),mean(TotalErrors(:,1)),'r*','LineWidth',1.5);
axis([-0.1 0.1 -0.1 0.1]);     xlabel('Cross-track');ylabel('Along-track');
grid on;

% 质量检验  消除孤立点
TotalErrors=TotalErrors(:,1:2);
Mean=mean(TotalErrors);
Delta=sqrt(sum((TotalErrors-Mean).^2,2));
Ind=find(Delta<2*mean(Delta));
TotalErrorsN=TotalErrors(Ind,:);
figure;scatter(TotalErrorsN(:,2),TotalErrorsN(:,1),'b.');
hold on;scatter(mean(TotalErrorsN(:,2)),mean(TotalErrorsN(:,1)),'r*','LineWidth',1.5);
xlabel('Cross-track');ylabel('Along-track');
axis([-0.1 0.1 -0.1 0.1]);
grid on;


