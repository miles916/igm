
clear all; clc; close all

%%%%%%%% Paths %%%%%%%
root='C:\Users\jouberto\Desktop\T&C';
path_func = [root '\TC_setups\Functions'];
addpath(path_func)
addpath(genpath([root '\TC_setups\Functions\topotoolbox'])) % Topotoolbox


SITE = 'Rolwaling'; % Catchment name
addpath(genpath([root, '/TC_Setups/',SITE,'/RUNS/INPUTS']));
simnm = 'IGM_test';

%%%% Load catchment grids %%%%%%%%%%%%%%%%

dtm_file = ['dtm_' SITE '_200m.mat']; % T&C pre-processing file
outlocation = [root,'/TC_outputs/',SITE,'/Distributed/',simnm,'/'];
mkdir(outlocation);
load(dtm_file); % Load T&C pre-processing grids (DTM, GLA, GLH, x, y, ...)

path_out = [outlocation SITE '_geologys.nc']; % Where to store the .nc file
DTM_Bedrock = DTM_orig-GLH; 

%%%%%% Provide a fake SMB maps (as would be obtained from T&C) %%%%%%%

ELA = 5600; % m a.s.l.
Accum_grad = 0.003;
Accum_abl = 0.015;

SMB = DTM.*0;
SMB(DTM > ELA) = Accum_grad.*(DTM(DTM > ELA) - ELA);
SMB(DTM <= ELA) = Accum_abl.*(DTM(DTM <= ELA) - ELA);
SMB(GLH == 0) = -10;

SMB_mean = nanmean(SMB(GLH > 0));

GLA_MAP2(MASK == 0) = 0;
GLH(MASK == 0) = 0;
DTM_Bedrock(MASK == 0) = 0;
DTM_orig(MASK == 0) = 0;

%%%% Create or update the .nc file %%%%%%%%%%%%

IGM_netcdf(path_out,x,y,flipud(GLA_MAP2),flipud(GLH),flipud(DTM_Bedrock), flipud(DTM_orig),flipud(SMB)) % Create .nc file

%%%% Run IGM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pyenv('Version', ... 
            'C:\Users\jouberto\.conda\envs\igm\python.exe', ... 
            'ExecutionMode','OutOfProcess')   % setup python virtual env

year_start = 2000;
year_end = 2020;
year_save = 1;  % How often to store IGM outputs (years)

igm_run_path = ['C:\Users\jouberto\Desktop\T&C\TC_setups\Rolwaling\' ...
    'Preprocessing\igm-run-rolwaling.py']; % path of IGM .py script
igm_arg = [' --geology_file ' path_out ...
            ' --working_dir ' outlocation ...
            ' --tstart ' num2str(year_start)...
            ' --tend ' num2str(year_end) ... 
            ' --tsave ' num2str(year_save)];

pyrunfile([igm_run_path igm_arg ]) % Run IGM

%%% IGM outputs

time = ncread([outlocation 'ex.nc'],'time');
thk = ncread([outlocation 'ex.nc'],'thk');
topg = ncread([outlocation 'ex.nc'],'topg');


%%
 
fi4 = figure('Renderer', 'painters', 'Position', [156.3333 86.3333 723.3334 512.0000]);
tiledlayout(2,2,'TileSpacing','compact')
nexttile
imagesc(flipud(reshape(thk(:,:,1),size(thk,2),size(thk,1)))); hold on;
set(gca,'Color',[0.8 0.8 0.8])
cb = colorbar; clim([0 400]); ylabel(cb,'Initial thickness [m]','FontSize',11)
cm = colormap('jet'); colormap([1,1,1; cm]);
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
nexttile
imagesc(flipud(reshape(thk(:,:,end),size(thk,2),size(thk,1)))); hold on;
set(gca,'Color',[0.8 0.8 0.8])
cb = colorbar; clim([0 400]); ylabel(cb,'Final thickness [m]','FontSize',11)
cm = colormap('jet'); colormap([1,1,1; cm]);
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
ax3 = nexttile;
imagesc(flipud(SMB)); hold on;
set(gca,'Color',[0.8 0.8 0.8])
cb = colorbar; clim([-4 4]); ylabel(cb,'Ice thickness change [m]','FontSize',11)
colormap(ax3,flipud(redblue));
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
ax4 = nexttile;
imagesc(flipud(reshape(thk(:,:,end)-thk(:,:,1),size(thk,2),size(thk,1)))); hold on;
set(gca,'Color',[0.8 0.8 0.8])
cb = colorbar; clim([-4 4]); ylabel(cb,'Ice thickness change [m]','FontSize',11)
colormap(ax4,flipud(redblue));
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])

nansum(SMB,'all')
nansum(thk(:,:,end)-thk(:,:,1),'all')


%% Figure with TopoZeko

% v = VideoWriter([outlocation 'IGM_3D'],'MPEG-4');
% v.FrameRate = 20;
% open(v)

figure
for tt = 1:size(thk,3)
    Surface = flipud(DTM_orig);
    GLH_3d = reshape(thk(:,:,tt),size(thk,2),size(thk,1));

    TopoZeko(fliplr(Surface-GLH_3d),fliplr(Surface),'extra_dimension',fliplr(GLH_3d),'vertical_scaling',0.5,'caxis',[0 500],...
    'D4_colormap','jet','D4_colormap_flipud','off','title',['Ice thickness [mm w.e.] - '...
    num2str(time(tt))],'view_orientation',[176.622638177387 70.8504878048781])
%     frame = getframe(gcf);
%     writeVideo(v,frame);
pause(0.1)
end 
% close(v)

% [16.6143336791176 47.7607317073171] : Kyzylsu
% [176.622638177387 70.8504878048781] : Langtang