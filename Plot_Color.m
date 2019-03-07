
function Plot_Color(NewObserveData0,Latitude,Longitude,Inflected_latitudes1,Inflected_longitudes1,TruthLatitude,TruthLongitude,CentLat, CentLon)
%画亮温的颜色图

R=1;
fanwei.Lat_min=max(floor(CentLat)-R,-89.9999);
fanwei.Lat_max=min(ceil(CentLat)+R,89.9999);
fanwei.Lon_min=max(floor(CentLon)-R,-179.9999);
fanwei.Lon_max=min(ceil(CentLon)+R,179.9999);
[X,Y] = meshgrid(1:size(NewObserveData0,2),1:size(NewObserveData0,1));
%调节插值的多少
Step=1/(1*4);
[X1,Y1] = meshgrid(1:Step:size(NewObserveData0,2),1:Step:size(NewObserveData0,1));



NewObserveData0([find(Latitude<fanwei.Lat_min); find(Latitude>fanwei.Lat_max);...
    find(Longitude<fanwei.Lon_min); find(Longitude>fanwei.Lon_max)]) = NaN;
Latitude([find(Latitude<fanwei.Lat_min); find(Latitude>fanwei.Lat_max)]) = NaN;%NaN
Longitude([find(Longitude<fanwei.Lon_min); find(Longitude>fanwei.Lon_max)]) = NaN;

NewLatitude = interp2(X, Y, Latitude, X1, Y1,'linear');
NewLongitude = interp2(X, Y, Longitude, X1, Y1,'linear');
NewObserveData = interp2(X, Y, NewObserveData0, X1, Y1,'cubic');
NewLatitude(:,end-(1/Step):end)=[];
NewLongitude(:,end-(1/Step):end)=[];
NewObserveData(:,end-(1/Step):end)=[];

%开始画图

figure('color','white');
latlim = [fanwei.Lat_min fanwei.Lat_max];
lonlim = [fanwei.Lon_min fanwei.Lon_max];
% latlim = [46 56];
% lonlim = [130 150];
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim, ...
      'Frame','on', 'MeridianLabel','on','ParallelLabel','on'); 
axis off;
setm(gca,'MLabelLocation',1);%
setm(gca,'PLabelLocation',1);%
% set(gca,'XTick',fanwei.Lat_min:0.5:fanwei.Lat_max);
% set(gca,'XTick',fanwei.Lon_min:0.5:fanwei.Lon_max);
Threshold=5;% 10 previous
dlevels = floor(min(min(NewObserveData))/Threshold)*Threshold : Threshold: ceil(max(max(NewObserveData))/Threshold)*Threshold;
for k = 1:length(dlevels) - 1
   NewObserveData(find(NewObserveData>dlevels(k) & NewObserveData<=dlevels(k+1))) = k;
end
NewObserveData(find(NewObserveData==dlevels(1))) = 1;
cmap = colormap(jet(length(dlevels) - 1));
colormap(cmap);
caxis([0 length(dlevels)-1]);
cbar = colorbar;
set(cbar,'Ticks',0:1:length(dlevels)-1,'TickLabels',dlevels);

pcolorm(NewLatitude, NewLongitude, NewObserveData);
geoshow(TruthLatitude,TruthLongitude,'Color','black');
setm(gca, 'fontsize', 16);
set(gcf, 'outerposition', get(0, 'screensize'));
% hold on;
% geoshow(Inflected_latitudes1,Inflected_longitudes1,'DisplayType','point');
% geoshow(Inflected_latitudes1,Inflected_longitudes1);

end