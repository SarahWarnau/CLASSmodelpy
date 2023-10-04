# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:43:54 2023

@author: sarah
"""

import cdsapi

c = cdsapi.Client()

# c.retrieve(
#     'reanalysis-uerra-europe-pressure-levels',
#     {
#         'variable': 'temperature',
#         'pressure_level': [
#             '500', '600', '700',
#             '750', '800', '825',
#             '850', '875', '900',
#             '925', '950', '975',
#             '1000',
#         ],
#         'year': '2011',
#         'month': '07',
#         'day': [
#             '01', '02', '03',
#             '04', '05', '06',
#             '07', '08', '09',
#             '10', '11', '12',
#             '13', '14', '15',
#             '16', '17', '18',
#             '19', '20', '21',
#             '22', '23', '24',
#             '25', '26', '27',
#             '28', '29', '30',
#             '31',
#         ],
#         'time': '12:00',
#         'format': 'netcdf',
#     },
#     'download.nc')

# c.retrieve(
#     'reanalysis-era5-pressure-levels-monthly-means',
#     {
#         'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
#         'variable': [
#             'geopotential', 'relative_humidity', 'specific_humidity',
#             'temperature', 'u_component_of_wind', 'v_component_of_wind', 
#             'vertical_velocity',
#         ],
#         'pressure_level': [
#             '500', '550', '600',
#             '650', '700', '750',
#             '775', '800', '825',
#             '850', '875', '900',
#             '925', '950', '975',
#             '1000',
#         ],
#         'year': '2011',
#         'month': '04',
#         'time': '12:00',
#         'area': [
#             44, -11, 36,
#             4,
#         ],
#         'format': 'netcdf',
#     },
#     'ERA5_Spain_201104.nc')

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
        'variable': [
            '100m_u_component_of_wind', '100m_v_component_of_wind', '10m_u_component_of_wind',
            '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature',
            'cloud_base_height', 'forecast_logarithm_of_surface_roughness_for_heat', 'surface_latent_heat_flux',
            'surface_pressure', 'surface_sensible_heat_flux', 'surface_solar_radiation_downward_clear_sky',
            'surface_thermal_radiation_downward_clear_sky', 'total_column_cloud_liquid_water', 'total_precipitation',
            'boundary_layer_height',
        ],
        'year': '2011',
        'month': '04',
        'time': '12:00',
        'area': [
            44, -11, 36,
            4,
        ],
        'format': 'netcdf',
    },
    'ERA5_Spain_201104_surf.nc')

# c.retrieve(
#     'reanalysis-era5-single-levels-monthly-means',
#     {
#         'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
#         'variable': 'geopotential',
#         'year': '2011',
#         'month': '04',
#         'time': '12:00',
#         'area': [
#            44, -11, 36,
#            4,
#         ],
#         'format': 'netcdf',
#     },
#     'geopot.nc')


#%%
# import	gsw
import	numpy	as	np
import	xarray	as	xr
# import	cmocean
import	matplotlib.pyplot	as	plt
import	cartopy.crs	as	ccrs
import	cartopy.feature	as	cfeature
from matplotlib.pyplot import cm

#%%

ds	=	xr.open_dataset('ERA5_Spain_201104.nc')
dssurf = 	xr.open_dataset('ERA5_Spain_201104_surf.nc')
geopot = xr.open_dataset('geopot.nc')
# P = ds.level*100
cp = 1004
cv = 717
kappa =  (cp - cv)/cp 
# pottemp = ds.t*(101325/P)**kappa

#%%
def section(da,	which={'latitude':40}):
				"""	returns	time	mean"""
				return	da.sel(**which,	method='nearest').mean('time')
# def plot_section(which,	title=None):
#     f,	ax	=	plt.subplots(2,2,	figsize=(15,10)) 
#     if	title	is not	None:
#         f.suptitle(title) 
#         cbkw	=	dict(orientation='horizontal',	pad=0.1) 
#         DS = section(ds,	which)
#     if 'longitude' in	which:		
#         x	=	DS.latitude
#     elif 'latitude' in	which:
#         x	=	DS.longitude
#         #	temperature
#         # 	
    
#     im	=	ax[0,0].pcolormesh(x.values, DS.level.values, DS.t*(1013.25/DS.level)**kappa,	cmap='Spectral_r') 
#     plt.colorbar(im,	ax=ax[0,0],	label='potential temperature [K]', **cbkw) 
#     ax[0,0].invert_yaxis()
#     ax[0,0].plot(x.values, section(dssurf, which).sp/100, c='k')
				#	salinity
# 				im	=	ax[0,1].pcolormesh(x,	S.level, 1000*section(S.salt,	which))
# 				plt.colorbar(im,	ax=ax[0,1],	label='salinity		[g/kg]', **cbkw)
				#	density
# 				im	=	ax[1,0].pcolormesh(x,	sigma0.level,	section(sigma0,	which),	cmap='inferno')
# 				plt.colorbar(im,	ax=ax[1,0],	label='potential	density		[kg/m^3]', **cbkw)
				#	stratification
# 				im	=	ax[1,1].pcolormesh(x,	N2.level,	section(1/np.sqrt(N2),	which)/3600,	cmap='magma',	vmax=1)
# 				plt.colorbar(im,	ax=ax[1,1],	label=r'buoyancy	period	1/N		[h]', **cbkw)


#%%
# lat = 40
# lonslice=slice(-1.5,2)
lons = np.arange(-0.75,0.8, 0.1)
lats = (-1/1.5)*lons + 40
pathlen = len(lons)

coords = - np.sqrt((lons*85.18)**2 + (lats*111.2)**2) + np.sqrt((0*85.18)**2 + (40*111.2)**2)  

x = xr.DataArray(lons, dims="path")
x['path'] = coords
y = xr.DataArray(lats, dims="path")
y['path'] = coords


dspath = ds.interp(longitude=x, latitude=y).mean('time')
dssurfpath = dssurf.interp(longitude=x, latitude=y).mean('time')
geopotpath = geopot.interp(longitude=x, latitude=y).mean('time')

#%%

# DS = section(ds.sel({'level':slice(600,1000)}).sel({'longitude':lonslice}),	{'latitude':lat})
DS = dspath.sel({'level':slice(600,1000)})
# DSsurf = section(dssurf.sel({'longitude':lonslice}),	{'latitude':lat})
DSsurf = dssurfpath
pottemp = DS.t*(1013.25/DS.level)**kappa
dp = -25#pottemp.isel({'level':slice(0,8)} ).level.values - pottemp.isel({'level':slice(1,9)} ).level.values
nlevels = len(DS.level.values)
dtheta = pottemp.isel({'level':slice(0, nlevels-1)} ).values - pottemp.isel({'level':slice(1,nlevels)} ).values
dthetadp = dtheta/dp
gamma = dthetadp*-1*9.81

# DS['theta'] = (('level', 'longitude'), np.array(pottemp, dtype='float32'))
DS['theta'] = (('level', 'path'), np.array(pottemp, dtype='float32'))

DSsurf['theta'] = ('path', np.array((DSsurf.t2m*(101325/DSsurf.sp)**kappa), dtype='float32'))
R = 287
DSsurf['rho'] = ('path', np.array(DSsurf.sp/(DSsurf.t2m*R), dtype='float32'))




# gamma2 = np.append(gamma, np.array([np.empty(lonlen)]), axis=0)
# gamma2 = np.append(gamma, np.array([np.empty(pathlen)]), axis=1)

# DS['gamma'] = (('level', 'longitude'), np.array(gamma2, dtype='float32'))
# DS['gamma'] = (('level', 'path'), np.array(gamma2, dtype='float32'))

# DS['sp'] = (('longitude'), np.array(DSsurf.sp.values, dtype='float32'))
DS['sp'] = (('path'), np.array(DSsurf.sp.values, dtype='float32'))


# n = len(DS.longitude.values)
n = len(DS.path.values)

# DS
#%%
def plot_section(which,	title=None):
    f,	ax	=	plt.subplots(1,1,	figsize=(15,10), dpi=300) 
    if	title	is not	None:
        f.suptitle(title) 
        cbkw	=	dict(orientation='horizontal',	pad=0.1) 
        # DS = section(ds.sel({'level':slice(800,1000)}).sel({'longitude':lonslice}),	which)
    if 'longitude' in	which:		
        x	=	DS.latitude
        ax.set_xlabel("latitude (deg. N)")
    elif 'latitude' in	which:
        x	=	DS.longitude
        ax.set_xlabel("longitude (deg. E)")
        #	temperature
        # 	
    
    im	=	ax.pcolormesh(x.values, DS.level.values, DS.theta,	cmap='Spectral_r') 
    plt.colorbar(im,	ax=ax,	label='potential temperature [K]', **cbkw) 
    CS1 = ax.contour(x.values, DS.level.values, DS.t, colors='k', levels=10)
    ax.clabel(CS1, inline=True, fontsize=10)
    
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.t,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='temperature [K]', **cbkw) 
    # CS2 = ax.contour(x.values, DS.level.values, DS.theta, colors='k', levels=12)
    
    # ax.clabel(CS2, inline=True, fontsize=10)
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.gamma,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='potential temperature gradient [K/m]', **cbkw) 
    
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.q*1000,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='specific humidity (g/kg)', **cbkw) 
    # CS1 = ax.contour(x.values, DS.level.values, DS.r, colors='k', levels=10)
    # ax.clabel(CS1, inline=True, fontsize=10)

    ax.set_ylabel("pressure (hPa)")
    ax.fill_between(DS.longitude, DS.sp/100, np.ones(n)*1013.25, color='grey', zorder=50)
    ax.set_xlim(DS.longitude.min(),DS.longitude.max())
    ax.set_ylim(DS.level.min(), DS.level.max())
    ax.invert_yaxis()
    
# plot_section({'latitude':lat},	title=f'Transect along {lat}N')


def plot_section_path(title=None):
    f,	ax	=	plt.subplots(1,1,	figsize=(15,10), dpi=300) 
    if	title	is not	None:
        f.suptitle(title) 
        cbkw	=	dict(orientation='horizontal',	pad=0.1) 
        # DS = section(ds.sel({'level':slice(800,1000)}).sel({'longitude':lonslice}),	which)
    # if 'longitude' in	which:		
    #     x	=	DS.latitude
    #     ax.set_xlabel("latitude (deg. N)")
    # elif 'latitude' in	which:
    #     x	=	DS.longitude
    #     ax.set_xlabel("longitude (deg. E)")
    #     #	temperature
        # 	
    x = DS.path
    
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.theta,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='potential temperature [K]', **cbkw) 
    # CS1 = ax.contour(x.values, DS.level.values, DS.t, colors='k', levels=10)
    # ax.clabel(CS1, inline=True, fontsize=10)
    
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.theta,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='potential temperature [K]', **cbkw) 
    # CS1 = ax.contour(x.values, DS.level.values, DS.t, colors='k', levels=10)
    # ax.clabel(CS1, inline=True, fontsize=10)
    
    im	=	ax.pcolormesh(x.values, DS.level.values, DS.t,	cmap='Spectral_r') 
    plt.colorbar(im,	ax=ax,	label='temperature [K]', **cbkw) 
    CS2 = ax.contour(x.values, DS.level.values, DS.theta, colors='k', levels=12)
    ax.clabel(CS2, inline=True, fontsize=10)
    
    # ax.clabel(CS2, inline=True, fontsize=10)
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.gamma,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='potential temperature gradient [K/m]', **cbkw) 
    
    # im	=	ax.pcolormesh(x.values, DS.level.values, DS.q*1000,	cmap='Spectral_r') 
    # plt.colorbar(im,	ax=ax,	label='specific humidity (g/kg)', **cbkw) 
    # CS1 = ax.contour(x.values, DS.level.values, DS.r, colors='k', levels=10)
    # ax.clabel(CS1, inline=True, fontsize=10)

    ax.set_ylabel("pressure (hPa)")
    ax.set_xlabel("Distance from 40N 0E (km)")
    ax.fill_between(DS.path, DS.sp/100, np.ones(n)*1013.25, color='grey', zorder=50)
    ax.set_xlim(DS.path.min(),DS.path.max())
    ax.set_ylim(DS.level.min(), DS.level.max())
    ax.invert_yaxis()
    
plot_section_path(title=f'Transect')

#%%
fig, ax = plt.subplots(dpi=300)
# n = len(DS.longitude.values)
# colors=cm.nipy_spectral(np.linspace(0,1,n))
colors=cm.nipy_spectral(np.linspace(0,1,pathlen))
# for i, lon in enumerate(DS.longitude.values[::-1]):
for i, path in enumerate(DS.path.values[::-1]):
    if i%2 == 0:
        # dslon = DS.sel({'longitude':lon})
        dspath = DS.sel({'path':path})
        # ax.plot(dslon.theta, dslon.level, label=lon, c=colors[i])
        ax.plot(dspath.theta, dspath.level, label=int(path), c=colors[i])
        # ax.scatter(dslon.sel({'level':dslon.sp/100}, method='nearest').theta, dslon.sp/100, color=colors[i])
ax.fill_between(np.sort(DS.sel({'level':DS.sp/100}, method='nearest').theta), -np.sort(-DS.sp/100), np.ones(n)*1013, color='grey', zorder=50)
ax.invert_yaxis()
ax.legend()
ax.set_xlabel("potential temperature [K]")
ax.set_ylabel("pressure (hPa)")

#%%
fig, ax = plt.subplots(dpi=300)
# n = len(DS.longitude.values)
# colors=cm.nipy_spectral(np.linspace(0,1,n))
colors=cm.nipy_spectral(np.linspace(0,1,pathlen))

# for i, lon in enumerate(DS.longitude.values[::-1]):
for i, path in enumerate(DS.path.values[::-1]):
    if i%2 == 0:
        # dslon = DS.sel({'longitude':lon})
        dspath = DS.sel({'path':path})
        # ax.plot(dslon.q*1000, dslon.level, label=lon, c=colors[i])
        ax.plot(dspath.q*1000, dspath.level, label=int(path), c=colors[i])
        # ax.scatter(dslon.sel({'level':dslon.sp/100}, method='nearest').q*1000, dslon.sp/100, color=colors[i])
        ax.scatter(dspath.sel({'level':dspath.sp/100}, method='nearest').q*1000, dspath.sp/100, color=colors[i])
# ax.fill_between(np.sort(DS.sel({'level':DS.sp/100}, method='nearest').q*1000), (DS.sp/100), np.ones(n)*1013, color='grey', zorder=50)
ax.invert_yaxis()
ax.legend()
ax.set_xlabel('specific humidity (g/kg)')
ax.set_ylabel("pressure (hPa)")


#%%


# dssel = ds.sel({'level':slice(750,1025), 'longitude':lonslice})
fig,ax = plt.subplots()
# ax.streamplot(x=dssel.longitude.values, y=dssel.level.values,
#                u=dssel.u.sel({"latitude":lat}).mean('time').values,
#                v=dssel.w.sel({"latitude":lat}).mean('time').values*1e2, density=1, linewidth=None)
# quiver = ax.quiver(dssel.longitude.values, dssel.level.values,
#           dssel.u.sel({"latitude":lat}).mean('time').values,
          # dssel.w.sel({"latitude":lat}).mean('time').values)
quiver = ax.quiver(dspath.path.values, dspath.level.values, 
                   dspath.v.values, dspath.w.values)
# ax.quiverkey(quiver, X=0.3, Y=1.1, U=10,
#              label='Quiver key, length = 10', labelpos='E')

# ax.fill_between(dssel.longitude.values, 
#                 dssurf.sel({'latitude':lat, 'longitude':lonslice}).sp.values.squeeze()/100, 
#                 np.ones(len(dssel.longitude.values))*1013, color='grey', zorder=50)
ax.fill_between(dspath.path.values, 
                dssurfpath.sp.values.squeeze()/100, 
                np.ones(pathlen)*1013, color='grey', zorder=50)

ax.invert_yaxis()
# ax.set_title(lat)

#%%

fig = plt.figure(figsize=(10, 5), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# ax.set_extent([-11, 4, 44, 36,], crs=ccrs.PlateCarree())
ax.set_extent([-2,2,39,42], crs=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.RIVERS, color='blue')
# ax.hlines(xmin=-1.5, xmax=2, y=lat, color='r')
# ax.plot([-1.5,0],[41,40], color='r')
ax.plot(DS.longitude, DS.latitude, color='r')
# ax.streamplot(dssurf.longitude.values, dssurf.latitude.values,
#           dssurf.u100.mean('time').values,
#           dssurf.v100.mean('time').values, color='black')
im = ax.contourf(dssurf.longitude.values, dssurf.latitude.values,
              geopot.mean('time').z/9.81, cmap='terrain', levels=50)
plt.colorbar(im, ax=ax, label='height (m)')
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
fig.tight_layout()

#%%

fig = plt.figure(figsize=(10, 5), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# ax.set_extent([-11, 4, 44, 36,], crs=ccrs.PlateCarree())
ax.set_extent([-6,4,36,44], crs=ccrs.PlateCarree())
# ax.coastlines(edgecolor='black')
ax.add_feature(cfeature.BORDERS, linewidth=1, edgecolor='black')
ax.add_feature(cfeature.OCEAN, zorder=10)


ax.streamplot(dssurf.longitude.values, dssurf.latitude.values,
          dssurf.u100.mean('time').values,
          dssurf.v100.mean('time').values, color='black', zorder=20)
im = ax.contourf(dssurf.longitude.values, dssurf.latitude.values,
              geopot.mean('time').z/9.81, cmap='terrain', levels=40)
plt.colorbar(im, ax=ax, label='height (m)')
ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='white', alpha=0.5, linestyle='--')
fig.tight_layout()


#%%

from scipy import optimize

def test_func(x, a, b, c):
    return a * np.sin(b * x) + c

# params, params_covariance = optimize.curve_fit(test_func, DS.path.values, geopotpath.z.values/9.81,
                                               # p0=[1300, 0.05, -20])
params = -1300, 0.03, 0
# print(params)
a, b, c = params
# a,b = params
plt.plot(DS.path, geopotpath.z/9.81)
plt.plot(DS.path.values, test_func(DS.path.values, a, b, c))


#%%

u = DSsurf.u100
v = DSsurf.v100

# Calculate wind speed
wind_speed = np.sqrt(u**2 + v**2)

plt.plot(u.path, wind_speed)
plt.plot(np.arange(0,55)*-1, -3.1/55 * np.arange(0,55) + 3.3)
# Calculate wind direction in degrees
# wind_direction_degrees = np.arctan2(v, u) * (180.0 / np.pi)

# # Make sure wind direction is between 0 and 360 degrees
# wind_direction_degrees = (wind_direction_degrees + 360) % 360


#%%

def find_nearest(array, value):
    """
    Function to find the index of the value in an array nearest to given value
    By Gijs van Leeuwen https://github.com/gjvanleeuwen

    Parameters
    ----------
    array : 1d.array
        1d np array wih random values, floats or ints
    value : float or int
        value to look for in array

    Returns
    -------
    idx : int
        index for nearest value

    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#%%
fig, ax = plt.subplots(dpi=300)
# n = len(DS.longitude.values)
# colors=cm.nipy_spectral(np.linspace(0,1,n))
colors=cm.nipy_spectral(np.linspace(0,1,pathlen))

# run1input.alpha = 0.33
# run1input.wg = 0.255
# run1input.w2 = run1input.wg
run1input.sw_test = True  



rho=1

# for i, lon in enumerate(DS.longitude.values[::-1]):
for i, path in enumerate(DS.path.values[::-1]):
    if i%2 == 0:
        # dslon = DS.sel({'longitude':lon})
        dspath = DS.sel({'path':path})
        # ax.plot(dslon.theta, dslon.level, label=lon, c=colors[i])
        ax.plot(dspath.theta, dspath.level, label=int(path), c=colors[i])
        dspathsurf = DSsurf.sel({'path':path})
        ax.scatter(dspathsurf.theta-3, dspathsurf.sp/100, color=colors[i], zorder=100)
        ax.scatter(dspathsurf.theta-3, dspathsurf.sp/100  - dspathsurf.blh*dspathsurf.rho*9.81/100, color=colors[i], zorder=100)
        ax.vlines( dspathsurf.theta-3, dspathsurf.sp/100 - dspathsurf.blh*dspathsurf.rho*9.81/100, dspathsurf.sp/100, color=colors[i], linestyle=':', zorder=100)
        # ax.scatter(dslon.sel({'level':dslon.sp/100}, method='nearest').theta, dslon.sp/100, color=colors[i])
        if path <= 0:
            
            x = path*-1000
            run1input.runtime = x/(run1input.U - (3.1/55000)*x)
            run1input.tstart = (12*3600 - run1input.runtime)/3600
            r1 = model(run1input)
            r1.run()
            x_i = find_nearest(r1.out.x, dspath.path.values*-1000)
            print(r1.out.x[x_i])
            ax.vlines(r1.out.theta[x_i], r1.out.P_h[x_i]/100, r1.out.Ps[x_i]/100, color=colors[i], zorder=100)
            
ax.fill_between(np.sort(DS.sel({'level':DS.sp/100}, method='nearest').theta), -np.sort(-DS.sp/100), np.ones(n)*1013, color='grey', zorder=50)
ax.invert_yaxis()
ax.legend()
ax.set_xlabel("potential temperature [K]")
ax.set_ylabel("pressure (hPa)")

#%%
from metpy import calc as mpcalc
from metpy.units import units

#%%

from runmodelMC import run1input
from modelMC import *
n_runs = 10 
# run1input.theta = 294
run1input.q = 8.5e-3
run1input.h = 600
times = np.linspace(4,11.9,n_runs)
run1input.alpha = 0.3
run1input.wg = 0.21#0.23
run1input.w2 = run1input.wg
run1input.sw_test = True  

fig, ax = plt.subplots(2,2, dpi=300, sharex=True)

for i, t in enumerate(times):
    run1input.tstart = t
    r1 = model(run1input)
    r1.run()
    t12 = find_nearest(r1.out.t, 12)
    ax[0,0].scatter(-r1.out.x[t12]/1000, r1.out.H[t12], color='r')
    ax[0,0].scatter(-r1.out.x[t12]/1000, r1.out.LE[t12], color='b')
    ax[0,1].scatter(-r1.out.x[t12]/1000, r1.out.h[t12], color='k')
    ax[1,0].scatter(-r1.out.x[t12]/1000, r1.out.theta[t12], color='k')
    ax[1,1].scatter(-r1.out.x[t12]/1000, r1.out.q[t12]*1000, color='k')
    ax[0,0].scatter(-r1.out.x[t12]/1000, r1.out.LE[t12] + r1.out.H[t12], color='k')
    
ax[0,0].plot(DSsurf.path, (DSsurf.sshf + DSsurf.slhf)/-3600, c='k', label='Q at 12h00')
ax[0,0].plot(DSsurf.path, DSsurf.slhf/-3600, c='b', label='LE at 12h00')
ax[0,0].plot(DSsurf.path, DSsurf.sshf/-3600, c='r', label='H at 12h00')
ax[0,1].plot(DSsurf.path, DSsurf.blh, c='r', label='H at 12h00')
ax[1,0].plot(DSsurf.path, DSsurf.theta, c='r', label='H at 12h00')
q = mpcalc.specific_humidity_from_dewpoint(pressure=DSsurf.sp*units.Pa, dewpoint=DSsurf.d2m*units.K)
ax[1,1].plot(DSsurf.path, q*1000, c='r', label='H at 12h00')

ax[0,0].set_ylabel('Q, H, LE')
# ax[0,1].set_ylabel('H')
ax[1,1].set_ylabel('q*1000')
ax[1,0].set_ylabel(r'$\theta$')
# ax[2,0].set_ylabel('RH(z=h)')
# ax[2,1].set_ylabel(r'T$_{2m}$')
# ax[0,2].set_ylabel('Q')
ax[0,1].set_ylabel('ABL top')
# ax[2,2].set_ylabel(r'$\Delta \theta$')
ax[0,0].set_xlabel('x')
    
fig.tight_layout()
# ax[0,0].legend()

#%%
