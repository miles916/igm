function IGM_netcdf_ESM(path_out,x,y,GLA,GLH,DEM_bedrock, DEM, U, V)

%%%%%%%% INPUT %%%%%%%%%%%%%%%%%
% NOTE - arrays need to be shaped such that x is dim1 and y is dim2; y
% coords should increase. ie from Matlab image coordinates, the arrays need
% to be flipped vertically, then transposed
% path_out = absolute path of where should be stored the .nc file
% x : vector of X coordinates in UTM
% y : vector of Y coordinates in UTM
% GLA: glacier mask array
% GLH: ice thickness array
% DEM_bedrock: bedrock DEM
% DEM: Surface DEM
% SMB : surface map balance map

%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%
% Create a the "geology.nc" file needed as input for IGM

row = size(GLH,2);
col = size(GLH,1);

% Create .nc file or open it if already existing
if exist(path_out, 'file') == 0
    ncid = netcdf.create(path_out,'NETCDF4');
else 
    delete(path_out);
    ncid = netcdf.create(path_out,'NETCDF4');
end 
try
    % Define dimensions
    dimid1 =  netcdf.defDim(ncid,'x',col);
    dimid2 =  netcdf.defDim(ncid,'y',row);
    % Define variables
    icemask_id = netcdf.defVar(ncid,"icemaskobs","NC_DOUBLE",[dimid1 dimid2]);
    topg_id = netcdf.defVar(ncid,"topg","NC_DOUBLE",[dimid1 dimid2]);
    usurf_id = netcdf.defVar(ncid,"usurfobs","NC_DOUBLE",[dimid1 dimid2]);
    thk_id = netcdf.defVar(ncid,"thkobs","NC_DOUBLE",[dimid1 dimid2]);
    uvel_id = netcdf.defVar(ncid,"uvelsurfobs","NC_DOUBLE",[dimid1 dimid2]);
    vvel_id = netcdf.defVar(ncid,"vvelsurfobs","NC_DOUBLE",[dimid1 dimid2]);
    
    x_id = netcdf.defVar(ncid,"x","NC_DOUBLE",dimid1);
    y_id = netcdf.defVar(ncid,"y","NC_DOUBLE",dimid2);
    % Add values to variables
    netcdf.putVar(ncid,icemask_id,GLA)
    netcdf.putVar(ncid,topg_id,DEM_bedrock)
    netcdf.putVar(ncid,usurf_id,DEM)
    netcdf.putVar(ncid,thk_id,GLH)
    netcdf.putVar(ncid,uvel_id,U)
    netcdf.putVar(ncid,vvel_id,V)
    netcdf.putVar(ncid,x_id,x)
    netcdf.putVar(ncid,y_id,y)
    % Close .nc file
    netcdf.close(ncid)
catch % if a failure, still need to close the file!
    % Close .nc file
    netcdf.close(ncid)
end