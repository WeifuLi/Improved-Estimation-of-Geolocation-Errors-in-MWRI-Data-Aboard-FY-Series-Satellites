function hdf_create(data)
%% Readdata

name=data.name;
Point_row_indices=data.Point_row_indices;
Choosen_columns=data.Choosen_columns;
Closest_latitudes=data.Closest_latitudes;
Closest_longitudes=data.Closest_longitudes;
Inflected_latitudes=data.Inflected_latitudes;
Inflected_longitudes=data.Inflected_longitudes;



%% Create a new file using the default properties.
NewName = [name(1:end-4), '1.HDF'];
FILE=NewName;

file = H5F.create (FILE, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

%
% Create group creation property list and set it to allow creation
% of intermediate groups.
%
gcpl = H5P.create ('H5P_LINK_CREATE');
H5P.set_create_intermediate_group (gcpl, 1);

%
% Create the group /date.  Note that /date do not
% exist yet.  This call would cause an error if we did not use the
% previously created property list.

group = H5G.create (file, '/Data', gcpl,'H5P_DEFAULT', 'H5P_DEFAULT');


%%
%write data
DATASET1='Point_row_indices';        %行索引
DATASET2='Choosen_columns';          %列索引
DATASET3='Closest_latitudes';        %真值点经度
DATASET4='Closest_longitudes';       %真值点纬度
DATASET5='Inflected_latitudes';      %插值点经度
DATASET6='Inflected_longitudes';     %插值点纬度

H1=length(Point_row_indices);       
% dims1=[H1,L1];       wdata1= zeros(dims1);

%
% Initialize data.
%
% for i=1: H1
%     for j=1: L1
%         ii=i-1;
%         jj=j-1;
%         wdata1(i,j) =  ii / (jj + 0.5) + jj;
%     end
% end

% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
%
space1 = H5S.create_simple (1,fliplr( H1), []);
%
% Create the dataset and write the floating podata to it.  In
% this example we will save the data as 64 bit little endian IEEE
% floating ponumbers, regardless of the native type.  The HDF5
% library automatically converts between different floating point
% types.
%

dset1 = H5D.create (group, DATASET1, 'H5T_STD_I32LE', space1, 'H5P_DEFAULT');
H5D.write (dset1, 'H5T_NATIVE_INT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Point_row_indices);


H2=length(Choosen_columns);  
space2 = H5S.create_simple (1,fliplr(H2), []);
dset2 = H5D.create (group, DATASET2, 'H5T_IEEE_F64LE', space2, 'H5P_DEFAULT');
H5D.write (dset2, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Choosen_columns);

H3=length(Closest_latitudes);         
space3 = H5S.create_simple (1,fliplr( H3), []);
dset3 = H5D.create (group, DATASET3, 'H5T_IEEE_F64LE', space3, 'H5P_DEFAULT');
H5D.write (dset3, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Closest_latitudes);

H4=length(Closest_longitudes);  
space4 = H5S.create_simple (1,fliplr( H4), []);
dset4 = H5D.create (group, DATASET4, 'H5T_IEEE_F64LE', space4, 'H5P_DEFAULT');
H5D.write (dset4, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Closest_longitudes);

H5=length(Inflected_latitudes);
space5 = H5S.create_simple (1,fliplr( H5), []);
dset5 = H5D.create (group, DATASET5, 'H5T_IEEE_F64LE', space5, 'H5P_DEFAULT');
H5D.write (dset5, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Inflected_latitudes);

H6=length(Inflected_longitudes);  
space6 = H5S.create_simple (1,fliplr( H6), []);
dset6 = H5D.create (group, DATASET6, 'H5T_IEEE_F64LE', space6, 'H5P_DEFAULT');
H5D.write (dset6, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT',Inflected_longitudes);

%
% Print all the objects in the files to show that intermediate
% groups have been created.
%
disp ('Objects in the file:');
H5O.visit (file, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', @op_func, []);

%
% Close and release resources.
%
H5P.close (gcpl);
H5G.close (group);
H5F.close (file);


end

function [status opdataOut]=op_func(rootId,name,opdataIn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Operator function for H5O.visit.  This function prints the
%name and type of the object passed to it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ('/');               % Print root group in object path %

objId=H5O.open(rootId,name,'H5P_DEFAULT');
type=H5I.get_type(objId);
H5O.close(objId);

%
% Check if the current object is the root group, and if not print
% the full path name and type.
%
if (~(name == '.'))         % Root group, do not print '.' %
    
    switch (type)
        case H5ML.get_constant_value('H5I_GROUP')
            disp ([name ' (Group)']);
        case H5ML.get_constant_value('H5I_DATASET')
            disp ([name ' (Dataset)']);
        case H5ML.get_constant_value('H5I_DATATYPE')
            disp ([name ' (Datatype)']);
        otherwise
            disp ([name ' (Unknown)']);
    end
    
    
end
status=0;
opdataOut=opdataIn;
end