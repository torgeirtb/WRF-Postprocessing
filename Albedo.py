# -*- coding: utf-8 -*-
"""
    Script for plotting albedo over Skibotn from wrfout files.
    
        Basemap features:
        - 'merc': conserves correct directions. Sylindrica projection.
            wiki: https://no.wikipedia.org/wiki/Mercators_projeksjon
        - llcrnrlon: longitude of lower left hand corner of the desired map domain (degrees).
        - llcrnrlat: latitude of lower left hand corner of the desired map domain (degrees).
        - urcrnrlon: longitude of upper right hand corner of the desired map domain (degrees).
        - urcrnrlat: latitude of upper right hand corner of the desired map domain (degrees).
        - resolutions: c - crude, l - low, i - intermediate, h - high, f - full
        
    Rieppi coordinates:
        - longitude: 20.6804  
        - latitude: 69.1867

 26.05.2016
"""


from netCDF4 import Dataset as NetCDFFile
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


wrf_path = '/Users/torgeir/Documents/wrfout'
geo_path = '/Users/torgeir/Documents/geo_em.d03.nc'

nc_wrf = NetCDFFile(wrf_path, mode='r')
nc_geo = NetCDFFile(geo_path, mode='r')

font = 15       # <--- Fontsize

# Extracting 3-D variables (time, south-north, west-east) from WRF
time = 2    # time of interest
lon = nc_wrf.variables['XLONG'][time,:,:]
lat = nc_wrf.variables['XLAT'][time,:,:]
albedo = nc_wrf.variables['ALBEDO'][time,:,:]


# Extracting top. height from geo-file.
lon_geo = nc_geo.variables['XLONG_M'][0,:,:]
lat_geo = nc_geo.variables['XLAT_M'][0,:,:]
hgt = nc_geo.variables['HGT_M'][0,:,:]


# Create a map over the Rieppi region (Lyngen alps)
map = Basemap(projection='tmerc', lon_0=20.6804,lat_0=69.1867,
              llcrnrlon=20,llcrnrlat=68.95,
              urcrnrlon=21.25,urcrnrlat=69.45,
              resolution='i')
map.drawcoastlines()        
map.drawmeridians(np.arange(19,22,0.5), labels=[0,0,0,1], fontsize=font)     
map.drawparallels(np.arange(69,70,0.25), labels=[1,0,0,0], fontsize=font)  


# Transform the lat/lon data to map coordinates
x,y = map(lon,lat)

# Wind barbs
xx = np.arange(0,len(x))
yy = np.arange(0,len(y))
points = np.meshgrid(xx,yy)

# Topography
xg, yg = map(lon_geo, lat_geo)
th = map.contour(xg, yg, hgt, 6, colors='k')  # Equipotential lines from 'HGT'

# Albedo as a contour-surface plot
cs = map.contourf(x, y, albedo ,cmap='greys')

# Rieppi coordinates
Rlon, Rlat = map(20.6804, 69.1867)

# ------ FIGURE --------
map.plot(Rlon, Rlat, 'ro', markersize=8)
cb = map.colorbar(cs, 'right', size='5%', pad='5%')
cb.ax.tick_params(labelsize=font) 
cb.set_label('Temp at 2m [$^{\circ}$ C]', rotation=270, fontsize=font)
cb.ax.tick_params(labelsize=font) 
map.drawmapscale(20.3, 69., 21.5, 69., 20, 
                 barstyle='fancy', fontsize=font, units='km')
plt.title('Domain 1, 01.12.14', fontsize=font)

# plt.savefig('../Figures/10mWindField_date.png', dpi=200, bbox_inches='tight')


