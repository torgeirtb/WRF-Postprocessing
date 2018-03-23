% ------------------------------------------------------------------------------------------
% Extracting netcdf variables from WRF-runs and saves variables of interest
% to a .mat-file.
% Explanation of the individual WRF variables can be found by entering
% "ncdisp(path)" in MATLAB. Or by "ncdump -h <filename>" in a terminal.
%
% NetCDF variables (subscripts denote domain numbers):
%   (u,v,w) : 3-D wind velocity
%   PHB : base-state geopotential
%   PH  : perturbation geopotential
%   T2  : Temperature at 2 meters a.g.l.
%   lon : degrees east
%   lat : degrees north
%   PSFC: Surface pressure [Pa]
%   HGT : Terrain height [m]
%   znw : height levels (eta-levels)
%   znu : half height levels
%
% Last edited: 21.March.2018, Torgeir
% ------------------------------------------------------------------------------------------

% Identify WRF-outfile and where to store extracted data
wrfpath = '/global/work/blasterdalen/WRFsimMasterThesis/WRF_simulations/WRF_Dec/';
yy = 2014;
mm = 12
dd = 21;
HH = 00;
MM = 00;
SS = 00;
domain = 3;
fileID   = strcat(wrfpath, 'wrfout_d0', num2str(domain), '_', ...
                  num2str(yy), '-', num2str(mm), '-', num2str(dd), '_', ...
                  num2str(HH), ':', num2str(MM), ':', num2str(SS));

outDataPath = '/global/work/blasterdalen/WRF_narvik/WRFnygfj/DataExtracts/'

% ----------- WRF variables ----------------------------------------------------------------
% 4-D variables: (west-east, north-south, height-level, time)
u  = ncread(fileID, 'U');
v  = ncread(fileID, 'V');
w  = ncread(fileID, 'W');
PHB= ncread(fileID, 'PHB');
PH = ncread(fileID, 'PH');

% 3-D variables: (west-east, south-north, time)
T2 = ncread(fileID, 'T2');
lon  = ncread(fileID, 'XLONG');
lat  = ncread(fileID, 'XLAT');
PSFC = ncread(fileID, 'PSFC');
HGT  = ncread(fileID, 'HGT');

% 2-D variables: (height-level, time)
znw  = ncread(fileID, 'ZNW');
znu  = ncread(fileID, 'ZNU');
% ------------------------------------------------------------------------------------------

% Simulation time;
time1          = ncread(fileID, 'Times')';
Timestring1    = datenum(time1)';

% Actual longitude and latitude of site
R_lon = 20.6804;
R_lat = 69.1867;

% Obtaining the latitude and longitude indexes of the Rieppi site
lonxy = lon(:, :, 1);
latxy = lat(:, :, 1);
[nx, ny] = size(squeeze(lonxy));

dist = zeros(nx, ny);
r = 6.370e6;        % Radius of the Earth
for i = 1:nx
   for j = 1:ny
        dist(i,j) = r*sqrt(cos((R_lat + latxy(i, j))/2)^2*(R_lon - lonxy(i, j))^2 ...
                    + (R_lat - latxy(i, j))^2);
   end
end

% Finding the indexes of the closest grid points
[~, index_min] = min(dist(:));
[lat_index, lon_index] = ind2sub(size(dist), index_min);

% Saving to .mat file
uu  = squeeze(u(lon_index, lat_index, :, :));
vv  = squeeze(v(lon_index, lat_index, :, :));
ww  = squeeze(w(lon_index, lat_index, :, :));
ph  = squeeze(PH(lon_index, lat_index, :, :));
phb = squeeze(PHB(lon_index, lat_index, :, :));
psfc= squeeze(PSFC(lon_index, lat_index, :));
t2  = squeeze(T2(lon_index, lat_index, :))';      %'
hgt = squeeze(HGT(lon_index1,lat_index1,:));

save(strcat(outDataPath, 'D0', num2str(domain), '_', num2str(dd), '_', ...
     num2str(mm), '.mat'), ...
     'Timestring', 'uu', 'vv', 'ww', 'ph', 'phb', 'psfc', 't2', 'znw', ...
     'znu', 'hgt', 'lonxy', 'latxy')

quit
