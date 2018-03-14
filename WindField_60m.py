# -*- coding: utf-8 -*-
"""
    Script for plotting wind field over Skibotn from wrfout files at 60 meters a.g.l.
    
        Basemap features:
        - 'tmerc': conserves correct directions. Sylindrica projection.
            wiki: https://no.wikipedia.org/wiki/Mercators_projeksjon
        - llcrnrlon: longitude of lower left hand corner of the desired map domain (degrees).
        - llcrnrlat: latitude of lower left hand corner of the desired map domain (degrees).
        - urcrnrlon: longitude of upper right hand corner of the desired map domain (degrees).
        - urcrnrlat: latitude of upper right hand corner of the desired map domain (degrees).
        - resolutions: c - crude, l - low, i - intermediate, h - high, f - full
        
    Rieppi coordinates:
        - longitude: 20.6804  
        - latitude: 69.1867

 11.05.2016
"""


from netCDF4 import Dataset as NetCDFFile
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# Fontsize 
font = 15

wrf_path = '/global/work/blasterdalen/WRF_Dec/wrfout_d03_2014-12-10_00:00:00'
geo_path = '/home/blasterdalen/Ymse/geo_em.d03.nc'

nc_wrf = NetCDFFile(wrf_path, mode='r')
nc_geo = NetCDFFile(geo_path, mode='r')

# Extracting 3-D variables (time, south-north, west-east) from WRF
time = 72    # time of interest
lon = nc_wrf.variables['XLONG'][time,:,:]
lat = nc_wrf.variables['XLAT'][time,:,:]
# 4-D variables (time, eta, lon, lat)
eta = 4 	# index of 4 gives (normally) height 60 meter a.g.l.
u   = nc_wrf.variables['U'][time,eta,0:300,0:300]
v   = nc_wrf.variables['V'][time,eta,0:300,0:300]

# Extracting top. height from geo-file
lon_geo = nc_geo.variables['XLONG_M'][0,:,:]
lat_geo = nc_geo.variables['XLAT_M'][0,:,:]
hgt = nc_geo.variables['HGT_M'][0,:,:]


# Create a map over the Rieppi region (Lyngen alps)
map = Basemap(projection='tmerc', lon_0=20.6804,lat_0=69.1867,
              llcrnrlon=20,llcrnrlat=68.94,
              urcrnrlon=21.25,urcrnrlat=69.45,
              resolution='h')
map.drawcoastlines()
map.drawmeridians(np.arange(19,22,0.5), labels=[0,0,0,1], fontsize=font)
map.drawparallels(np.arange(69,70,0.25), labels=[0,0,0,0], fontsize=font)


# Transform the lat/lon data to map coordinates
x,y = map(lon,lat)

# Wind barbs
xx = np.arange(0,len(x), 2) 	# Numpy arange (start, stop, step)
yy = np.arange(0,len(y), 2)
points = np.meshgrid(xx,yy)

# Topography
xg, yg = map(lon_geo, lat_geo)
th = map.contour(xg, yg, hgt, 6, colors='k')  # Equipotential lines from 'HGT'


WindField = np.sqrt(u*u + v*v)
cs = map.contourf(x, y, WindField,cmap='summer')
map.barbs(x[points], y[points], u[points], v[points], barbcolor='r')

# Rieppi coordinates
Rlon, Rlat = map(20.6804, 69.1867)

# ------ FIGURE ---------
map.plot(Rlon, Rlat, 'bo', markersize=8)
cb = map.colorbar(cs, 'right', size='5%', pad='5%')
cb.ax.tick_params(labelsize=font)
map.drawmapscale(20.3, 68.98, 20.6804, 69.1867, 20, 
                 barstyle='fancy', fontsize=font)
cb.set_label('Wind speed [m/s]', rotation=270, fontsize=font)
plt.title('Domain 3, 10.12.14, 12:00', fontsize=font)
plt.savefig('/home/blasterdalen/Figures/WindField_u1005_10dec.png', dpi=200, bbox_inches='tight')




