function [HDF,Errors,PixelErrors]=ImprovedGER(DATA,RegionName,Method)
% Latitude ����
% Latitude γ��
% ObserveData �����¶�
% LandSeaMask ��½����
% TruthLatitude TruthLongitude GSHHS�����ĸ߾��ȵĺ�����
Latitude=DATA.Latitude;
Longitude=DATA.Longitude;
ObserveData=DATA.ObserveData;
LandSeaMask=DATA.LandSeaMask;
TruthCoastline=DATA.TruthCoastline;

TruthLatitude=TruthCoastline(:,1);
TruthLongitude=TruthCoastline(:,2);
TruthIndex=TruthCoastline(:,3);


% ��Ҫ���������
Filename=DATA.Filename;
Point_row_indices=[];
Choosen_columns=[];
Inflected_latitudes=[];
Inflected_longitudes=[];
Closest_latitudes=[];
Closest_longitudes=[];
Errors=[];
PixelErrors=[];

addpath(genpath('libicp'));
ResultPath=[RegionName '���ӻ����_' Filename '\'];
mkdir(ResultPath);
%Step 1: ��landseamask������ɸѡ��½�ͺ���߽�ĺ����ߵ�
[T5,T6]=ImagePreprocessing(LandSeaMask);
% figure;imshow(T5,[]);
% saveas(gcf,[ResultPath, Method 'T5.tif']);
% figure;imshow(T6,[]);
% saveas(gcf,[ResultPath, Method 'T6.tif']);
%�������γ����
%T5=T5(:,301:end-300);
T5=T5';
%T6=T6(:,301:end-300);
T6=T6';
%ObserveData=ObserveData(:,301:end-300);
ObserveData=ObserveData';
%Latitude=Latitude(:,301:end-300);
Latitude=Latitude';
%Longitude=Longitude(:,301:end-300);
Longitude=Longitude';

% ���н�
LatCenter=(Latitude(:,127)+Latitude(:,128))/2;
LonCenter=(Longitude(:,127)+Longitude(:,128))/2;
TrackDir=[diff(LatCenter) diff(LonCenter)];
HeadingAngle=acos(abs(TrackDir(:,1))./sqrt(TrackDir(:,1).^2+TrackDir(:,2).^2));
HeadingAngle=[HeadingAngle;HeadingAngle(end)];

%Step 2: ѡȡT6��ÿһ������ֵ��Ӧ��FOV����뾶R�����FOV���������������ߵ�����
R=10;
Miss=0;
for i=200+1:1.5*R:size(T6,1)-200
    for j=R+2:1.5*R:size(T6,2)-R-1
        
       %% Coastline points dtection
        Ix1=i-R/2:i+R/2;
        Iy1=j-R/2:j+R/2;
        T=T6(Ix1,Iy1);
        CentLat=Latitude(i,j);
        CentLon=Longitude(i,j);        
        if sum(sum(T))<=1
            continue;
        end
        Data=ObserveData(Ix1,Iy1);
        % Step 3: �жϰ뾶ΪR/2���������±仯���ֵ�Ƿ����15
        M0=[-0.5,0,0.5];C0 = conv2(Data,M0,'valid');C0=padarray(C0,[0 1]);[Max0,Ti0]=max(abs(C0(:)/9));
        M1=[0 0 -0.5;0 0 0; 0.5 0 0];C1 = conv2(Data,M1,'valid');C1=padarray(C1,[1 1]);[Max1,Ti1]=max(abs(C1(:)/sqrt(306)));
        M2=[0.5;0;-0.5];C2 = conv2(Data,M2,'valid');C2=padarray(C2,[1 0]);[Max2,Ti2]=max(abs(C2(:)/15));
        M3=[-0.5 0 0;0 0 0;0 0 0.5];C3 = conv2(Data,M3,'valid');C3=padarray(C3,[1 1]);[Max3,Ti3]=max(abs(C3(:)/sqrt(306)));
        [MaxALL,Indmax]=max([Max0,Max1,Max2,Max3]);
        if MaxALL<0.8
            continue;
        end
        Ti=[Ti0,Ti1,Ti2,Ti3];Size=[size(C0);size(C1);size(C2);size(C3)];
        Ti=Ti(Indmax);Size=Size(Indmax,:);[Tix,Tiy]=ind2sub(Size,Ti);
        T1=padarray(T,[1 1]);
        Tem=T1(Tix:Tix+2,Tiy:Tiy+2);
        if sum(Tem(:))==0
            continue;
        end
        %         figure;imshow(T6,[]);hold on; plot(j,i,'ro','MarkerSize',10);
        Ix1=i-R:i+R;
        Iy1=j-R:j+R;
        Lat=Latitude(Ix1,Iy1);
        Lon=Longitude(Ix1,Iy1);
        Data1=ObserveData(Ix1,Iy1);
        Coastline=T5(Ix1,Iy1);
        %         figure;plot3(Lat,Lon,Data1,'.','MarkerSize',10);xlabel('Latitude'); ylabel('Longitude');zlabel('Brightness Temperature');
        % Step 4: ��MWRI���ݵĵ���ֱ���9 km*15 km���в�ֵ
        Res=1;
        Step1=1/(Res*15);Step2=1/(Res*9);
        [X,Y] = meshgrid(1:size(Data1,2),1:size(Data1,1));
        [X1,Y1] = meshgrid(1:Step2:size(Data1,2),1:Step1:size(Data1,1));
        Coastline= interp2(X,Y,double(Coastline),X1,Y1,'nearest');
        %���Բ�ֵ���ȡ�γ��
        Lat1 = interp2(X,Y,Lat,X1,Y1);
        Lon1 = interp2(X,Y,Lon,X1,Y1);
        %������������¶�
        Data_cubic = interp2(X,Y,Data1,X1,Y1,'cubic');
        g1=SurfaceFitting(Data_cubic,Lat1,Data1,Step1,Step2);
        %         figure;plot3(Lat1,Lon1,g1,'.','MarkerSize',10);  xlabel('Latitude'); ylabel('Longitude');zlabel('Brightness Temperature');
        % Step 5: �ҵ���ֵ��������¶ȱ仯������
        tic
        [Dx,Dy]=gradient(g1);
        D=sqrt(Dx.^2+Dy.^2);
        T=1.2/Res;
        D2=D>T;
        D2=bwareaopen(D2,200*Res*Res);
        D2=bwmorph(D2,'open',5*Res);
        [Bw,~]=bwlabel(D2);
        D3=D2&Coastline; BwValue=unique(Bw(D3>0)); BwValue=setdiff(BwValue,0);
        if isempty(BwValue)==1
            continue;
        end
        BwRegionValue=BwChosenValue(Bw,BwValue);
        
        BwRegionValue=imfill(BwRegionValue,'holes');%         figure;imshow(BwRegionValue,[]);
        PadSize=10;
        BwRegionValue= padarray(BwRegionValue,[PadSize PadSize],'replicate');
        BwRegionValue=bwmorph(BwRegionValue,'thin',Inf);
        BwRegionValue=BwRegionValue(PadSize+1:end-PadSize,PadSize+1:end-PadSize);
        % ȥë��
        BwRegionValue=BwRemoveSpur(BwRegionValue);
        %         figure;imshow(BwRegionValue,[]);
        BwRegionValue=bwareaopen(BwRegionValue,60*Res,8);
        toc
        if sum(sum(BwRegionValue))<60*Res
            continue;
        end
        %         BwRegionValue=BwMaxArea(BwRegionValue);
        % ��ÿ����Ĵ��߷������ҵ����±仯���ĵ�
        [Result]=MaxValueinVertical(BwRegionValue,D2,D);
        
        
        %% Geolocation error measurement
        %         figure;imshow(Result,[]);
        [IxR,IyR]=find(Result==1);
        Centroid=round([mean(IxR), mean(IyR)]);
        DistC=pdist2(Centroid,[IxR,IyR]);
        MajorAxisLength=max(DistC);
        
        %%    Ѱ�Ҹ�������ֵ����������
        Latx=Lat1(Centroid(1),Centroid(2));
        Lonx=Lon1(Centroid(1),Centroid(2));
        Latx1=Lat1(Centroid(1)+1,Centroid(2)+1);
        Lonx1=Lon1(Centroid(1)+1,Centroid(2)+1);
        DisR=MajorAxisLength*sqrt((Latx1-Latx).^2+(Lonx1-Lonx).^2);
        %         Dis=zeros(length(TruthLatitude),1);
        Error=(repmat([Latx Lonx],length(TruthLatitude),1)-[TruthLatitude TruthLongitude]);
        Error=sqrt(Error(:,1).^2+Error(:,2).^2);
        [minErr,Indmin]=min(Error);
        IndT=find(TruthIndex==TruthIndex(Indmin));
        Ind=find(Error<DisR);
        Ind=intersect(IndT,Ind);
        if isempty(Ind)==1
            continue;
        end
        %         Error1=Error(Ind)
        SLat=TruthLatitude(Ind);
        SLon=TruthLongitude(Ind);
        %         figure; hold on;plot(SLat,SLon,'r.');hold on; plot(Lat1(IndD2),Lon1(IndD2),'g.');
        %         figure;worldmap world; geoshow(SLat,SLon);hold on; geoshow(Lat1(IndD2),Lon1(IndD2),'DisplayType','point');
        p=[Lat1(Result==1) Lon1(Result==1)];% Detected coastline
        q=[SLat SLon];%Actual coastline
        if isempty(q)==0
            Dist2=pdist2(p,q);
            Theshold=max(10*sqrt((Latx1-Latx).^2+(Lonx1-Lonx).^2),5*minErr);
            Dist21=Dist2<=Theshold;
            %             Dist21=BwMaxArea(Dist21);
            Dist21=bwareaopen(Dist21,500);
            % �ҵ�����ÿ�����ĺ����ߵ��������ʵ�����ߵ�
            Dist22=Dist2;
            Dist22(~Dist21)=1;
            [~,MinInd]=min(Dist22,[],2);
            IndMin=sub2ind(size(Dist2),[1:size(Dist2,1)]',MinInd);
            MinRegion=zeros(size(Dist2));
            MinRegion(IndMin)=1;
            Se=strel('line',20,0);
            MinRegion=imdilate(MinRegion,Se);
            MinRegion=MinRegion&Dist21;
            MinRegion=bwmorph(MinRegion,'dilate',8);
            [DistIndx,DistIndy]=find(MinRegion==1);
            %             [DistIndx,DistIndy]=find(Dist21==1);
            DistIndx=unique(DistIndx);
            DistIndy=unique(DistIndy);
            if length(DistIndy)<10
                continue;
            end
            p=p(DistIndx,:);
            q=q(DistIndy,:);
            % ICP�㷨�ҵ����ĺ����ߺ���ʵ�ĺ����ߵ�ƽ�Ʊ任
            MinPQ=min([q; p]);
            q=q-MinPQ;
            p=p-MinPQ;
            p=1*p';q=1*q';
            %             [Qx,~]=gradient(q);
            %             Q=sqrt(Qx(1,:).^2+Qx(2,:).^2);
            %             BW=find(Q<0.05);
            %             BW=BwMaxArea(BW);
            %             IndB=find(BW>0);
            %             q=q(:,IndB);
            Tr_fit = icpMex(q,p,eye(3),-1,'point_to_plane');
            if Tr_fit(1,1)>0.98&&abs(Tr_fit(1,3))<0.2&&abs(Tr_fit(2,3))<0.2
                % ���ӻ���ʾ�㼯��׼���
                T_fit  = Tr_fit(1:2,1:2)*p + Tr_fit(1:2,3)*ones(1,size(p,2));
                Dist2=pdist2(T_fit',q');
                [~,MinInd]=min(Dist2,[],2);
                IndMin=sub2ind(size(Dist2),[1:size(Dist2,1)]',MinInd);
                MinRegion=zeros(size(Dist2));
                MinRegion(IndMin)=1;
                MinRegion=MinRegion&(Dist2<0.01);
                [DistIndx,DistIndy]=find(MinRegion==1);
                if length(DistIndx)<0.5*size(p,2)
                    continue;
                end
                T_fit=T_fit(:,DistIndx);
                q=q(:,DistIndy);
                
%                 figure,axis equal,hold on; ms=8; lw=2; fs=16;
%                 plot(p(2,:)+MinPQ(2),p(1,:)+MinPQ(1),'b','MarkerSize',ms,'LineWidth',lw);
%                 plot(q(2,:)+MinPQ(2),q(1,:)+MinPQ(1),'r','MarkerSize',ms,'LineWidth',lw);
%                 legend('Detected costline', 'True costline');
% %                 plot(T_fit(2,:)+MinPQ(2),T_fit(1,:)+MinPQ(1),'.b','MarkerSize',ms,'LineWidth',lw);
%                 set(gca,'FontSize',fs);
                p=p(:,DistIndx);
                Error=q-p;
                Error=mean(Error,2);
                p=p+MinPQ';               
                q=q+MinPQ';
                Cnt=0;Cnt1=0;Cnt2=0;
%                 figure;imshow(Result,[]);
                for Row=1:1/Step1:size(Result,1)
                    Cnt=Cnt+1;
                    Point_row_indice=Ix1(Cnt);
                    TemRow=Result(Row,:);
                    if sum(TemRow)==0
                        continue;
                    else
                        Ind_y=find(TemRow==1);
                        Ind_y=Ind_y(1);
                        Cnt1=Cnt1+1;
                    end
%                     hold on; plot(Ind_y,Row,'ro','MarkerSize',ms);
                    Choosen_column=Iy1(1)-1+Ind_y*Step2;
                    Inflected_latitude=Lat1(Row,Ind_y);
                    Inflected_longitude=Lon1(Row,Ind_y);
                    [IsE,Ind]=ismember([Inflected_latitude,Inflected_longitude],p','rows');
                    if IsE==0
                        continue;
                    else
                        Cnt2=Cnt2+1;
                        Closest_latitude=q(1,Ind);
                        Closest_longitude=q(2,Ind);
%                         hold on; plot(Inflected_longitude,Inflected_latitude,'go','MarkerSize',ms,'LineWidth',lw);
%                         hold on; plot(Closest_longitude,Closest_latitude,'gx','MarkerSize',ms,'LineWidth',lw);
                    end
                    
                    
                    Point_row_indices=[Point_row_indices Point_row_indice];
                    Choosen_columns=[Choosen_columns Choosen_column];
                    X=Point_row_indice;Y=fix(Choosen_column);
                    Along_dir=norm([Latitude(X+1,Y)-Latitude(X,Y) Longitude(X+1,Y)-Longitude(X,Y)]);
                    Cross_dir=norm([Latitude(X,Y+1)-Latitude(X,Y) Longitude(X,Y+1)-Longitude(X,Y)]);                    
                    Inflected_latitudes=[Inflected_latitudes Inflected_latitude];
                    Inflected_longitudes=[Inflected_longitudes Inflected_longitude];
                    %                 Closest_latitudes=[Closest_latitudes Inflected_latitude+Tr_fit(1,3)];
                    %                 Closest_longitudes=[Closest_longitudes Inflected_longitude+Tr_fit(2,3)];
                    Closest_latitudes=[Closest_latitudes Closest_latitude];
                    Closest_longitudes=[Closest_longitudes Closest_longitude];
                    CrossTrack=(Inflected_longitude-Closest_longitude)*sin(HeadingAngle(Point_row_indice))+(Inflected_latitude-Closest_latitude)*cos(HeadingAngle(Point_row_indice));
                    AlongTrack=-(Inflected_longitude-Closest_longitude)*cos(HeadingAngle(Point_row_indice))+(Inflected_latitude-Closest_latitude)*sin(HeadingAngle(Point_row_indice));
%                     Errors=[Errors;  Inflected_longitude-Closest_longitude Inflected_latitude-Closest_latitude];
                    Errors=[Errors;  CrossTrack AlongTrack];
                    
                   PixelErrors=[PixelErrors; CrossTrack/Cross_dir AlongTrack/Along_dir];
                    
                end                            
%                 hold on; plot(p(2,1:8:end),p(1,1:8:end),'bo','MarkerSize',ms,'LineWidth',lw);
%                 hold on; plot(q(2,1:8:end),q(1,1:8:end),'ro','MarkerSize',ms,'LineWidth',lw);
%                 legend('Detected costline', 'True costline');
%                 xlabel('Longitude'); ylabel('Latitude');
%                 title(['ѡȡ��ĸ���' num2str(Cnt2)  ';  ƽ��: lat�� ' num2str(Error(1),2) '; lon��' num2str(Error(2),2)]);
%                 saveas(gcf,[ResultPath, num2str(i), '_',num2str(j) '��׼ͼ' Method '.tif'  ]);
                          
                % ���ӻ���ʾ�����¶Ƚ��
%                 Plot_Color(ObserveData,Latitude,Longitude,Lat1(Result==1),Lon1(Result==1),TruthLatitude,TruthLongitude, CentLat, CentLon);
%                 saveas(gcf,[ResultPath, num2str(i), '_',num2str(j) '����ͼ' Method '.tif'  ]);
            else
                Miss=Miss+1;
            end
        end
        close all;
    end
end
% ����HDF�ĸ�ʽ�����ڶ�λ����
Miss
Errors=Errors(:,1:2);
Mean=mean(Errors);
Delta=sqrt(sum((Errors-Mean).^2,2));
Ind=find(Delta<2*mean(Delta));
Errors=Errors(Ind,:);
PixelErrors=PixelErrors(Ind,:);
HDF.name=Filename;
HDF.Point_row_indices=int32(Point_row_indices(Ind));           %�к�
HDF.Choosen_columns=Choosen_columns(Ind);            %�к�
HDF.Closest_latitudes=Closest_latitudes(Ind);        %����㾭��
HDF.Closest_longitudes=Closest_longitudes(Ind);      %�����γ��
HDF.Inflected_latitudes=Inflected_latitudes(Ind);              %��ֵ�㾭��
HDF.Inflected_longitudes=Inflected_longitudes(Ind);            %��ֵ��γ��
end
