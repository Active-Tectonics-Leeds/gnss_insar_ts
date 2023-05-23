#!/usr/bin/python
## P.Espin 
### PLOT time series and change the file of GPS to LOS

# import libraries
import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib.transforms as mtransforms
import math 

#Directory path
directorio= os.getcwd()

## Files
dfDesc_2016 = pd.read_csv(directorio+'/Desc_2016_time_new.csv')
dfAsc_2016 = pd.read_csv(directorio+'/Asc_2016_time_new.csv')
## GPS files
dfGPS = pd.read_csv(directorio+'/VC1.csv')

## Change to dates
dfDesc_2016['Dates'] = pd.to_datetime(dfDesc_2016['Dates'])
dfGPS['Dates'] = pd.to_datetime(dfGPS['Dates'])
dfAsc_2016['Dates'] = pd.to_datetime(dfAsc_2016['Dates'])

### Change to LOS
## Frame 142D_09148
inc=39.5621
head=-167.92738
sininc=math.sin(inc)
cosinc=math.cos(inc)
sinhead=math.sin(head)
coshead=math.cos(head)

###


fig=plt.figure(figsize=(20,15))
ax = fig.add_subplot(111, polar=True)

ax1 = plt.subplot(3,1,1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax1.text(0.0, 0.1, "a)", transform=ax1.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
plt.plot(dfDesc_2016.Dates, dfDesc_2016.N, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
ax1.set_title("Descending", fontsize=16, fontweight='bold')
plt.legend(loc='upper left', fontsize='small');
ax1.xaxis.set_visible(False)


ax2 = plt.subplot(3,1,2, sharex=ax1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax2.text(0.0, 0.1, "b)", transform=ax2.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
plt.plot(dfAsc_2016.Dates, dfAsc_2016.N, color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
ax2.set_title("Ascending", fontsize=16, fontweight='bold')
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.ylabel('LOS Displacement (mm)', fontsize=18, fontweight='bold', x= -5)  
ax1.xaxis.set_visible(False)

ax3 = plt.subplot(3,1,3, sharex=ax1)
trans = mtransforms.ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
ax3.text(0.0, 0.1, "c)", transform=ax3.transAxes + trans,fontsize='large', verticalalignment='top',     bbox=dict(facecolor='white', edgecolor='none', pad=3.0))
### Plot the GPS in LOS check the equation 
plt.plot(dfGPS.Dates, ((((dfGPS.dN*sinhead)-(dfGPS.dE*coshead))*sininc)+(dfGPS.dU*cosinc)), color='blue', marker="o", label='N', linestyle='None', markersize=6, linewidth=0.5)
ax3.set_title("GPS (LOS)", fontsize=16, fontweight='bold')
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.xlabel('Dates', fontsize=18,fontweight='bold')  
plt.ylabel('LOS Displacement (mm)', fontsize=18, fontweight='bold', x= -5)  
fig.suptitle("Time Series", fontweight='bold', fontsize=18, y=0.98)


plt.tight_layout()
plt.savefig('Time_Series_LOS.jpg', format='jpg', dpi=400, bbox_inches='tight')
plt.show()


