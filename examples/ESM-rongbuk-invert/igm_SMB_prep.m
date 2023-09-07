%master_ContinuitySMB_ex.m - Master script to estimate surface mass balance
%distribution for a glacier from inputs of ice thickness, thinning, and
%velocity, based on the continuity equation (see, e.g. Bisset et al, 2020: https://doi.org/10.3390/rs12101563)
%
% New glaciers require this template to be coded with paths for the full set of
% inputs, and any additional needed preprocessing. 
%
% Other m-files required: C2xyz.m, index_nanfill.m, remove_small_zones.m,
% segment_Gmask_slope2.m, subset_geo.m, through_fluxes.m, zonal_aggregate.m
%
% Author: Evan Miles
% Work address: Swiss Federal Research Institute WSL
% Email: evan.miles@wsl.ch
% Sep 2019; Last revision: 03-June-2020

%% ------------- BEGIN CODE --------------
clear 
close all

%set home directory
homedir = 'C:\Users\miles\Documents\GitHub\igm\examples\ESM-rongbuk-invert';
srcpath = 'C:\Users\miles\Documents\GitHub\grid_continuity_SMB\data\Rongbuk';
addpath(genpath('C:\Users\miles\Documents\GitHub\grid_continuity_SMB\code'))
addpath('C:\Users\miles\Documents\GitHub\igm\')
addpath(homedir)
%glacier ID
P.Glacier = 'Rongbuk'

%title of output directory
datatitle = ['test_' P.Glacier];

%inputs        
V.pathx = fullfile(srcpath,'HMA_G0120_vx.tif');
V.pathy = fullfile(srcpath,'HMA_G0120_vy.tif');
V.mult=1; %scale to convert input units to m/a
DH.path = fullfile(srcpath,'15.09991_dH.tif');
DH.mult=1; %scale to convert input units to m/a
THX.path = fullfile(srcpath,'RGI60-15.09991_thickness_composite.tif');
DEM.path = fullfile(srcpath,'15.09991_AW3D.tif');

%settings
P.plotouts=1; %output plots or not
P.exports=1; %save geotiffs or not
P.DX = 100; %resolution to run calculations at 
P.segdist=300; %effective linear distance between flowbands
P.Vfilter=0; %swtich to smooth velocity data or not
P.dhfilter=0; %switch to peform 2x 3-sigma outlier removal (from overall dataset - only if erroneous pixels are common)
P.THXfilter=0; %switch to apply simply Gaussian filter to thickness data (useful for field thickness measurements)
P.umult=0; %switch for column-average velocity [0.8-1] is physical range. '0' estimates based on THX dist. '2' estimates for each pixel. 
P.Vreproj=0; %0 if velocities are provided oriented in the correct coordinate system, 1 if they are in the source data projection, [2 if they are true north], [3 to determine from slope]
N.fdivfilt=2; %use VanTricht gradient filters (2), just flux filter (1) or not at al (0)
N.uncertainty=0; %0 for 'simple run without uncertainty, otherwise N runs

% initialization
addpath(genpath([homedir '\code'])) %add path to related scripts

outpath = fullfile(homedir,[num2str(P.DX) 'mgrid_' date]);
mkdir(outpath)
    
    %% Load all input data, resample to common grid, etc
    [THX,DEM,V,DH]=load_subset_all(P,THX,DEM,V,DH);
    
    %% Cleaning raw data if needed
    
    %identify likely errors
    V.iERR=(abs(V.Uraw)>400)|(abs(V.Vraw)>400);
    if P.Vfilter==1 %gaussian low-pass filter removing extreme variations
        V.Uraw=imgaussfilt(V.Uraw,5);
        V.Vraw=imgaussfilt(V.Vraw,5);
    end

    if P.THXfilter==1% filter thickness data
        THX.dataR=THX.data;
        THX.data=imgaussfilt(THX.dataR,4); %gaussian low-pass filter. Important for thickness maps derived from field data
        THX.data(MASK0==0)=0;
    end
    
    if P.dhfilter==1 %if noisy or gappy, apply some preprocessing to the dH
        DH.data2=DH.data;
        DH.errthresh1=3.*nanstd(DH.data(:));
        DH.data2(abs(DH.data2)>DH.errthresh1)=NaN;
        DH.errthresh2=3.*nanstd(DH.data2(:));
        DH.data2(abs(DH.data2)>DH.errthresh2)=NaN;
        DH.dH3=imgaussfilt(DH.data2); %smooths dH slightly, expands NaN around bad DH data
    else
        DH.dH3=DH.data;
    end
    
    %% REPROJECT AND RESAMPLE ALL INPUTS (to THX coordinate system), derive MASK, etc
    N = resample_inputs(P.DX,THX,V,DEM,DH);
    N.MASK=N.THX>0;

    N.DEM(isnan(N.DEM))=nanmin(N.DEM(:));
    %% fill gaps in reprojected inputs if needed
    
    %gap-fill dH based on elevation
    N.DH0=N.DH; %unfilled values
    N.DH= index_nanfill(N.DH,N.DEM);
    N.DH((N.MASK==0))=0;
    
    N.U(isnan(N.U))=0;
    N.V(isnan(N.V))=0;

% %% create ncfile for igm - DEPRECATED
% % s=ncinfo('observation.nc')
% 
% nccreate('observation.nc','y',"Dimensions",{"y",length(N.y3)})
% ncwrite("observation.nc",'y',flipud(N.y3'))
% nccreate('observation.nc','x',"Dimensions",{"x",length(N.x3)})
% ncwrite("observation.nc",'x',N.x3)
% 
% nccreate('observation.nc','usurfobs',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'usurfobs',flipud(N.DEM)')
% nccreate('observation.nc','thkobs',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'thkobs',flipud(N.THX)')
% nccreate('observation.nc','icemaskobs',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'icemaskobs',double(flipud(N.MASK))')
% 
% nccreate('observation.nc','uvelsurfobs',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'uvelsurfobs',double(flipud(N.U))')
% nccreate('observation.nc','vvelsurfobs',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'vvelsurfobs',double(flipud(N.V))')
% 
% nccreate('observation.nc','thkinit',"Dimensions",{"x","y"})
% ncwrite("observation.nc",'thkinit',flipud(N.THX)')
% 
% ncdisp('observation.nc')
[x,y]=pixcenters(N.Rout,size(N.DEM));
y=fliplr(y);%flip everything out of image coords - note that TF runs in [x,y] coordinates - flip and transpose necessary for conversion from image referencing
cd(outpath)
IGM_netcdf_ESM(fullfile(outpath,'observation3.nc'),x,y,flipud(int8(N.MASK))',flipud(N.THX)',flipud(N.DEM-N.THX)',flipud(N.DEM)',flipud(N.U)',flipud(N.V)'); %flips are to switch from image space to cartesian ordering

%% IGM post
ncf=fullfile(homedir,'geology-optimized.nc');
nci=ncinfo(ncf);
thx_opt=flipud(ncread(ncf,'thk')');
geotiffwrite('THX_opt.tif',thx_opt,N.Rout,'GeoKeyDirectoryTag',THX.info.GeoTIFFTags.GeoKeyDirectoryTag,'TiffTags' , struct ( 'Compression' , Tiff.Compression.Deflate))