#!/usr/bin/python
## P.Espin 
### PLOT time series and change the file of GPS to LOS
### Jess Payne update 19/05/23:
# Equation to calculate LoS GNSS displacements from three component GNSS displacements (Fialko & Simons, 2001;
# Stephens et al, 2020)
# The below script takes a .tenv3 three component GNSS file as downloaded from Nevada Geodetic GPS Portal
# and convert the three components into LoS direction to compare to LoS displacments calculated using
# Sentinel-1 InSAR.
# Linear fits to GNSS and InSAR data are calculated and plotted.
# .tenv3 parameter extraction and linear vertical GNSS fit translated from John Elliott MATLAB script
#
# InSAR data is calculated using LiCS processing tools.
#
# Inputs required:
# 1. .h5/.nc file as output from LiCSBAS or licsar2licsbas
# 1b. Parameter file as output from LiCSBAS or licsar2licsbas (noramlly output as EQA.dem_par)
# 2. .tenv3 file downloaded from http://geodesy.unr.edu/NGLStationPages/gpsnetmap/GPSNetMap.html
# 3. Lat, lon and name of GNSS site (find on website in 2.)
# 4. InSAR LiCS frame name

# import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.transforms as mtransforms
import math
import netCDF4
from scipy.stats import linregress
import h5py
from datetime import datetime

# this is Andrew Watson's library of functions, see https://github.com/Active-Tectonics-Leeds/interseismic_practical
import interseis_lib as lib

#%%
#Directory path
directorio= '/nfs/a285/homes/eejap/plots/gnss_insar_ts/'

filename = './data/MXTM.NA.tenv3'
# nc_file = '/nfs/a285/homes/eejap/reunwrap_tests/078A_07049_131313/078A_07049_131313_gold_casc_ml1_clip.nc'
h5_file = './data/078A_07049_131313_gucmT.ml1.clip.h5'
par_file_asc = './data/078A_07049_131313_gucmT_ml1_clip.dem_par'
formatSpec = '%4C%8s%10f%6f%5f%2f%7f%7f%10f%10f%10f%6f%10f%8f%9f%9f%9f%10f%10f%f%[^\n\r]'

gnss_lon = -98.9790
gnss_lat = 19.4840
buffer_radius = 10
insar_frame = "078A_07049_131313"
gnss_station = "MXTM"

dfgnss = pd.read_csv(filename, delimiter=r"\s+")
# nc_dataset = netCDF4.Dataset(nc_file)
# #%% View parameters of InSAR data in netCDF file
# print('Dimensions:')
# for dim in nc_dataset.dimensions:
#     print(dim, len(nc_dataset.dimensions[dim]))
# #%%
# print('Variables:')
# for var in nc_dataset.variables:
#     print(var, nc_dataset.variables[var].shape)

#%% get imdates
with h5py.File(h5_file, 'r') as file:
    imdates = file['imdates']
    imdates = imdates[:]  
    vel_asc = file['vel']
    vel_asc = vel_asc[:]
    cum_asc = file['cum']
    cum_asc = cum_asc[:]

#%% complete using h5 file
# read array dimensions from par file
width_asc = int(lib.get_par(par_file_asc,'width'))
length_asc = int(lib.get_par(par_file_asc,'nlines'))

# width_desc = int(lib.get_par(par_file_desc,'width'))
# length_desc = int(lib.get_par(par_file_desc,'nlines'))

# get corner positions
corner_lat_asc = float(lib.get_par(par_file_asc, 'corner_lat'))
corner_lon_asc = float(lib.get_par(par_file_asc,'corner_lon'))

# corner_lat_desc = float(lib.get_par(par_file_desc,'corner_lat'))
# corner_lon_desc = float(lib.get_par(par_file_desc,'corner_lon'))

# get post spacing (distance between velocity measurements)
post_lat_asc = float(lib.get_par(par_file_asc,'post_lat'))
post_lon_asc = float(lib.get_par(par_file_asc,'post_lon'))

# post_lat_desc = float(lib.get_par(par_file_desc,'post_lat'))
# post_lon_desc = float(lib.get_par(par_file_desc,'post_lon'))

# calculate grid spacings
lat_asc = corner_lat_asc + post_lat_asc*np.arange(1,length_asc+1) - post_lat_asc/2
lon_asc = corner_lon_asc + post_lon_asc*np.arange(1,width_asc+1) - post_lon_asc/2

# lat_desc = corner_lat_desc + post_lat_desc*np.arange(1,length_desc+1) - post_lat_desc/2
# lon_desc = corner_lon_desc + post_lon_desc*np.arange(1,width_desc+1) - post_lon_desc/2



#%% convert imdates to good format
dates = []
for date_num in imdates:
        date_str = str(date_num)      
        date_obj = datetime.strptime(date_str, "%Y%m%d")
        dates.append(date_obj)

#%% Read lat and lon from nc file

#lat = nc_dataset.variables['lat'][:]
#lon = nc_dataset.variables['lon'][:]
#cum = np.flip(nc_dataset.variables['cum'], axis=0)
#vel = np.flip(nc_dataset.variables['vel'][:,:], axis=0)



#%%

# Define the extent of the image using latitude and longitude values
lat_min, lat_max = lat_asc.min(), lat_asc.max()
lon_min, lon_max = lon_asc.min(), lon_asc.max()
# Plot the 'vel' variable with latitude and longitude on the axes
plt.imshow(vel_asc, extent=[lon_min, lon_max, lat_min, lat_max])

# Add a marker for the specific point (gnss_lat, gnss_lon)
plt.plot(gnss_lon, gnss_lat, 'rx', markersize=8)

# Add labels and title to the plot
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('{}'.format(insar_frame))
colorbar = plt.colorbar()
colorbar.set_label('LoS Velocity (mm/yr)')
plt.savefig('./outputs/insar_los_map.jpg', dpi=400, bbox_inches='tight')
plt.show()

#%% Find the indices of the nearest grid cell to the poi
lat_index = np.abs(lat_asc - gnss_lat).argmin()
lon_index = np.abs(lon_asc - gnss_lon).argmin()

#%%
# Extract the subset of data within the buffer
# cum_ts = np.flip(cum_asc[:, lat_index, lon_index], axis=0) # if using nc file
cum_ts = cum_asc[:, lat_index, lon_index] # if using cum file

# Extract the time dimension
# time = nc_dataset.variables['time'][:]
# #%%
# # Plot the buffer data
# plt.plot(time, cum_ts)
# plt.xlabel('Time')
# plt.ylabel('Cumulative displacement (LoS, mm)')
# #plt.title('Variable Values within Buffer')
# plt.show()

#%%
# Define the desired date range
start_date = "20170615"
end_date = "20190504"

# Define the desired date range
start_date = datetime.strptime(start_date, "%Y%m%d")
end_date = datetime.strptime(end_date, "%Y%m%d")

# Find the indices of the 'dates' array that correspond to the desired date range
start_index = next(idx for idx, t in enumerate(dates) if t >= start_date)
end_index = next(idx for idx, t in enumerate(dates) if t <= end_date)
end_index = next((idx for idx, date in enumerate(dates) if date > end_date), len(dates))
end_index -= 1

#%% Convert InSAR dates to decimals for calculations
# Taken from https://github.com/sczesla/PyAstronomy/blob/master/src/pyasl/asl/decimalYear.py
dates_dec = []
for d in dates:
    year = d.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)
    yearElapsed = (d) - (startOfThisYear)
    yearDuration = (startOfNextYear) - (startOfThisYear)
    fraction = yearElapsed/yearDuration
    date_dec = year + fraction
    dates_dec.append(date_dec)
    
# Extract the subset of 'cum' data for the desired date range
cum_subset = cum_ts[start_index:end_index+1]
dates_subset = dates[start_index:end_index+1]
dates_dec_subset = dates_dec[start_index:end_index+1]

# Plot the subset of 'cum' data
plt.plot(dates_subset, cum_subset)

# Add labels and title to the plot
plt.xlabel('Time')
plt.ylabel('Cumulative')
plt.title('Cumulative Data')

#%% Reformat to decimal year, E, N, U, remove mean and rescale to mm
dfdata = pd.DataFrame({
    'Column2': dfgnss.iloc[:, 2],
    'Column8': dfgnss.iloc[:, 8],
    'Column10': dfgnss.iloc[:, 10],
    'Column12': dfgnss.iloc[:, 12],
    'Column14': dfgnss.iloc[:, 14],
    'Column15': dfgnss.iloc[:, 15],
    'Column16': dfgnss.iloc[:, 16]
})
#%%
dfdata_array = np.squeeze(np.array(dfdata))
mean_value = np.mean(dfdata_array, axis=0).squeeze().reshape(1,-1)
binary_array = np.array([0, 1, 1, 1, 0, 0, 0])
scaling_array = np.array([1, 1000, 1000, 1000, 1000, 1000, 1000])

mean_value = np.tile(mean_value, (dfdata_array.shape[0], 1))
binary_array = np.tile(binary_array, (dfdata_array.shape[0], 1))
scaling_array = np.tile(scaling_array, (dfdata_array.shape[0], 1))

dfdata_array = scaling_array * (dfdata_array - binary_array * mean_value)
ndays = len(dfdata_array)
date = dfdata_array[:, 0]
E = dfdata_array[:, 1]
N = dfdata_array[:, 2]
U = dfdata_array[:, 3]
eE = dfdata_array[:, 4]
eN = dfdata_array[:, 5]
eU = dfdata_array[:, 6]
wE = 1 / eE ** 2
wN = 1 / eN ** 2
wU = 1 / eU ** 2

column_names = ['Dates', 'dN', 'dE', 'dU', 'Sn', 'Se', 'Su']
dfGPS=pd.DataFrame(dfdata_array, columns=column_names)

#%%
## Change to dates
#dfGPS['Dates'] = pd.to_datetime(dfGPS['Dates'], format='%Y.%j')
dfGPS['Dates'] = pd.to_datetime(dfGPS['Dates'], format='%Y') + pd.to_timedelta((dfGPS['Dates'] % 1) * 365, unit='D')

#%%
### Change to LOS
## Frame 078A_07049_131313
## Take inc and heading from LiCS Portal metadata.txt for frame of interest
inc=39.6507*(np.pi/180)
head=-10.804014
az=(360+head)*(np.pi/180)
sininc=math.sin(inc)
cosinc=math.cos(inc)
sinaz=math.sin(az)
cosaz=math.cos(az)
GPS_dLOS = (((dfGPS.dN*sinaz)-(dfGPS.dE*cosaz))*sininc)+(dfGPS.dU*cosinc)

#%%

fig=plt.figure(figsize=(20,15))
ax = fig.add_subplot(111, polar=True)

ax1 = plt.subplot(3,1,1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax1.text(0.0, 0.1, "a)", transform=ax1.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
ax1.set_title("GNSS Vertical Displacement ({})".format(gnss_station), fontsize=16)
plt.plot(dfGPS.Dates, dfGPS.dU, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)

ax2 = plt.subplot(3,1,2, sharex=ax1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax2.text(0.0, 0.1, "b)", transform=ax2.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
### Plot the GPS in LOS check the equation 
plt.plot(dfGPS.Dates, GPS_dLOS, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
ax2.set_title("GNSS LoS Displacement ({})".format(gnss_station), fontsize=16)
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)


ax3 = plt.subplot(3,1,3, sharex=ax1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax3.text(0.0, 0.1, "c)", transform=ax3.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
plt.plot(dates_subset, cum_subset, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
ax3.set_title("Ascending ({}) InSAR LoS Displacement".format(insar_frame), fontsize=16)
#ax3.xaxis.set_visible(False)

# ax2 = plt.subplot(3,1,1)
# trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
# ax2.text(0.0, 0.1, "a)", transform=ax1.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
# plt.plot(dates, cum_ts, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
# ax2.set_title("Descending", fontsize=16, fontweight='bold')
# plt.legend(loc='upper left', fontsize='small');
# ax2.xaxis.set_visible(False)

plt.xlabel('Dates', fontsize=18,fontweight='bold') 
plt.ylabel('mm', fontsize=18, x= -5)   
fig.suptitle("Displacement time-series at point {}, {}".format(gnss_lon,gnss_lat), fontweight='bold', fontsize=18, y=0.98)

plt.tight_layout()
plt.savefig('./outputs/{}_{}_disp_ts.jpg'.format(gnss_station, insar_frame), format='jpg', dpi=400, bbox_inches='tight')
plt.show()

#%%
## Calculate linear fit for vertical GNSS
G = np.column_stack((date, np.ones(ndays)))

Q = np.concatenate((np.linalg.pinv(G.T.dot(np.diag(wE)).dot(G)),
                    np.linalg.pinv(G.T.dot(np.diag(wN)).dot(G)),
                    np.linalg.pinv(G.T.dot(np.diag(wU)).dot(G))), axis=1)

me = np.column_stack((Q[0:2, 0:2].dot(G.T.dot(np.diag(wE)).dot(E)),
                     Q[0:2, 2:4].dot(G.T.dot(np.diag(wN)).dot(N)),
                     Q[0:2, 4:6].dot(G.T.dot(np.diag(wU)).dot(U))))

me_plot= str(round(me[0,2], 2))
Q_plot = str(round(np.sqrt(Q[0,4]), 2))
me_line = me[0,2] * date + me[1,2]

## Calculate linear regression for LoS GNSS
slope, intercept, r_value, p_value, std_err = linregress(date, GPS_dLOS)
slope_plot = format(slope, '.2f')
std_err_plot = format(std_err, '.2f')
regression_line = slope * date + intercept

## Calculate linear regression for LoS InSAR
slope_i, intercept_i, r_value_i, p_value_i, std_err_i = linregress(dates_dec_subset, cum_subset)
slope_plot_i = format(slope_i, '.2f')
std_err_plot_i = format(std_err_i, '.2f')
regression_line_i = slope_i * date + intercept_i

#%% Plot time-series with linear-fit
fig=plt.figure(figsize=(20,15))
ax = fig.add_subplot(111, polar=True)

ax1 = plt.subplot(3,1,1)
#trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax1.text(0.0, 0.1, "a)", transform=ax1.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
ax1.set_title("GNSS Vertical Displacement ({})".format(gnss_station), fontsize=16)
plt.plot(dfGPS.Dates, dfGPS.dU, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
plt.plot(dfGPS.Dates, (me[0,2] * date + me[1,2]), 'g')
plt.text(dfGPS.Dates.iloc[400], 150, 'Linear fit: ' + me_plot + ' +/- ' + Q_plot + ' mm/yr', fontsize=16)

ax2 = plt.subplot(3,1,2, sharex=ax1)
#trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax2.text(0.0, 0.1, "b)", transform=ax2.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
### Plot the GPS in LOS check the equation 
plt.plot(dfGPS.Dates, GPS_dLOS, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
plt.plot(dfGPS.Dates, regression_line, 'g')
plt.text(dfGPS.Dates.iloc[400], 100, 'Linear fit: ' + slope_plot + ' +/- ' + std_err_plot + ' mm/yr', fontsize=16)
ax2.set_title("GNSS LoS Displacement ({})".format(gnss_station), fontsize=16)


ax3 = plt.subplot(3,1,3, sharex=ax1)
#trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax3.text(0.0, 0.1, "c)", transform=ax3.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
plt.plot(dates_subset, cum_subset, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
plt.plot(dfGPS.Dates, regression_line_i, 'g')
plt.text(dfGPS.Dates.iloc[400], -400, 'Linear fit: ' + slope_plot_i + ' +/- ' + std_err_plot_i + ' mm/yr', fontsize=16)
ax3.set_title("Ascending ({}) InSAR LoS Displacement".format(insar_frame), fontsize=16)

plt.xlabel('Dates', fontsize=18,fontweight='bold') 
plt.ylabel('mm', fontsize=18, x= -5)   
fig.suptitle("Displacement time-series at point {}, {}".format(gnss_lon,gnss_lat), fontweight='bold', fontsize=18, y=0.98)

plt.tight_layout()
plt.savefig('./outputs/{}_{}_disp_ts_lin_fit.jpg'.format(gnss_station, insar_frame), dpi=400, bbox_inches='tight')
plt.show()

#%% Start all time-series at 0 to compare
GPS_dLOS_first = GPS_dLOS.iloc[0]
GPS_dLOS_zero = GPS_dLOS - GPS_dLOS_first

dfGPS_dU_first = dfGPS.dU.iloc[0]
dfGPS_dU_zero = dfGPS.dU - dfGPS_dU_first

cum_first = cum_subset[0]
cum_zero = cum_subset - cum_first

#%% Calculate new linear fit for zeroed vels
## Zero linear fits
me_line_zero = me_line - me_line[0]
regression_line_zero = regression_line - regression_line[0]
regression_line_i_zero = regression_line_i - regression_line_i[0]

#%% plot 'normalised' time-series
plt.figure(figsize=(10,6))
plt.plot(dfGPS.Dates, dfGPS_dU_zero, label="Vertical GNSS ({})".format(gnss_station), marker="o", linestyle='None', color='blue', markersize=3)
plt.plot(dfGPS.Dates, GPS_dLOS_zero, label="LoS GNSS ({})".format(gnss_station), marker="o", linestyle='None', color='cornflowerblue', markersize=3)
plt.plot(dates_subset, cum_zero, label="Reunwrapped LoS InSAR ({})".format(insar_frame), marker="o", linestyle='None', color='green', markersize=3)

# Plot zeroed linear fits
plt.plot(dfGPS.Dates, me_line_zero, label="Vertical GNSS Linear Fit ({})".format(gnss_station), color = 'blue')
plt.plot(dfGPS.Dates, regression_line_zero, label="LoS GNSS Linear Fit ({})".format(gnss_station), color = 'cornflowerblue')
plt.plot(dfGPS.Dates, regression_line_i_zero, label="Reunwrapped LoS InSAR Linear Fit ({})".format(insar_frame), color = 'green')

# Plot linear velocity
plt.text(dfGPS.Dates.iloc[325], -20, 'Vertical GNSS velocity: ' + me_plot + ' +/- ' + Q_plot + ' mm/yr', fontsize=10)
plt.text(dfGPS.Dates.iloc[325], -40, 'LoS GNSS velocity: ' + slope_plot + ' +/- ' + std_err_plot + ' mm/yr', fontsize=10)
plt.text(dfGPS.Dates.iloc[325], -60, 'LoS InSAR velocity: ' + slope_plot_i + ' +/- ' + std_err_plot_i + ' mm/yr', fontsize=10)

plt.xlabel('Dates') 
plt.ylabel('mm') 
plt.title("Displacement time-series at point {}, {}".format(gnss_lon,gnss_lat), fontweight='bold')
plt.legend(loc='lower left')

plt.savefig('./outputs/{}_{}_disp_ts_lin_fit_zeroed.jpg'.format(gnss_station, insar_frame), dpi=400, bbox_inches='tight')


