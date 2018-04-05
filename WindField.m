%:-------------------------------------------------------------------------
% Testscript for visualizing simulations outputs from WRF for one domain. 
% 
% Explanation of the individual WRF variables can be found by entering
% "ncdisp(path)" in MATLAB. Or by "ncdump -h <filename>" in a terminal.
% NOTE: The order of arguments in variables might not be the same when
% using "ncdump" in the terminal or "ncdisp" in MATLAB. It is therefore a
% good advice to use "ncdisp" to see MATLAB's order of arguments for each
% variable. 
% 
% Some useful NetCDF variables from WRF output files:
%   XLAT  : degrees north, South is negative. 
%           (west-east, south-north, time). [degrees]
%   XLONG : degrees east, West is negative.             --
%   ZNU   : half height levels (i.e. half eta-levels). (bottom-top, time)
%   ZNW   : height levels (i.e. eta-levels).                   -- 
%   U     : x-wind component 
%           (west-east - staggered, south-north, eta, time). [m/s]
%   V     : y-wind component 
%           (west-east, south-north - staggered, eta, time). [m/s]
%   W     : z-wind component      --
%   PH    : perturbation geopotential (lon, lat, eta, time). [m2/s2]
%   PHB   : base-state geopotential             --
%   T2    : Temperature at 2 meters (lon, lat, time). [K]
%   PSFC  : Surface pressure (lon, lat, time). [Pa]
% 
% Some useful variables from geographical data:
%   XLAT_M: Latitude (west-east, south-north, time). [degrees]
%   XLONG_M: Longitude              --
%   HGT_M : Topography height (west-east, south-north, time). [meters MSL]
% 
% Last edited: 05.April.2018, Torgeir
%:-------------------------------------------------------------------------

close all 
clear all 
clc

addpath ../MatlabFunctions/m_map
addpath ../MatlabFunctions/gebco
addpath ../MatlabFunctions

eta = 4; 	% index of 4 gives (normally) height 60 meter a.g.l.
timestep = 72;

% Path and file name
wrfID = '../WRFoutputs/wrfout_d01_2014-01-01_12:00:00';
geoID = '../WRFoutputs/geo_em.d01.nc';

% Read WRF variables 
LAT = ncread(wrfID, 'XLAT');
LON = ncread(wrfID, 'XLONG');
U = ncread(wrfID, 'U');
V = ncread(wrfID, 'V');

% Reduce redundant dimensions
lat = squeeze(LAT(:, :, timestep));
lon = squeeze(LON(:, :, timestep));
u = squeeze(U(:, :, eta, timestep));
v = squeeze(V(:, :, eta, timestep));

% Interpolate to theta-points
u_theta = NaN(size(u, 1) - 1, size(u, 2));
v_theta = NaN(size(v, 1), size(v, 2) - 1);
for i = 1:size(u, 1) - 1
    u_theta(i, :) = 0.5*(u(i, :) + u(i + 1, :));
    v_theta(i, :) = 0.5*(v(:, i) + v(:, i + 1));
end
windMag = sqrt(u_theta.^2 + v_theta.^2);

% Read geographical variables (these should be semi-static and thus we can
% read from time index 1.
GEO_lat = ncread(geoID, 'XLAT_M');
GEO_lon = ncread(geoID, 'XLONG_M');
GEO_hgt = ncread(geoID, 'HGT_M');

geo_lat = squeeze(GEO_lat(:, :, 1));
geo_lon = squeeze(GEO_lon(:, :, 1));
geo_hgt = squeeze(GEO_hgt(:, :, 1));

% Tidy up 
clear GEO* i LAT LON u U v V 


% ---- Plotting ----
% Start/end of domain
min_lon = min(min(geo_lon));
max_lon = max(max(geo_lon));
min_lat = min(min(geo_lat));
max_lat = max(max(geo_lat));

fig1 = figure(1);
hold all

% Projection
m_proj('lambert', 'long', [double(min_lon) double(max_lon)], ...
       'lat', [double(min_lat) double(max_lat)]);

% Wind magnitude with colorbar
[~, l] = m_contourf(lon, lat, windMag); set(l, 'linestyle', 'none');
caxis([min(min(windMag)) max(max(windMag))])
cb = colorbar('v');
set(get(cb, 'ylabel'), 'string', 'Wind speed [m/s]', 'fontsize', 12)

% Plot wind quivers for visualizing directions
m_quiver(lon, lat, u_theta, v_theta, 'color', 'r')

% Point of interest
[xR, yR] = m_ll2xy(20.6804, 69.1867);
plot(xR, yR, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10)

% Coastline (bathymetry at 0 m from GEBCO)
load gebco0.mat gebco0
region = [-8 10 78 82];
gebco0(gebco0(:, 2) < min_lon | gebco0(:, 2) > max_lon | ...
       gebco0(:, 1) < min_lat | gebco0(:, 1) > max_lat) = NaN;
m_line(gebco0(:, 2), gebco0(:, 1), 'color', 'k', 'linewidth', 1.2)

% Topography
m_contour(geo_lon, geo_lat, geo_hgt, 0:100:2000, 'color', 'k', ...
          'ShowText','off', 'linewidth', .6)
m_contour(geo_lon, geo_lat, geo_hgt, 0:500:1500, 'color', 'k', ...
          'ShowText','on', 'linewidth', .9)

% Grid and lon/lat - annotations
m_grid('box', 'on', 'xtick', ...
       double(round(min_lon)):2:double(round(max_lon)), ...
       'tickdir', 'out', 'ytick', ...
       double(round(min_lat)):double(round(max_lat)), ...
       'linest', ':', 'tickstyle', 'dd')

% Window size
set(fig1, 'Position', [0 0 800 800])
% save2pdf('../Figures/TestfigFrom_PlottingOutputVariables.pdf')



