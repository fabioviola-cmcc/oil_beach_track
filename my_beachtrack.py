#!/usr/bin/python

# requirements
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import traceback
import numpy as np
import imageio
import pdb


# input/output configuration
inputFile = '/Users/fviola/code/medslik_netcdf_utils/fusion_track_beach/input/spill_properties.nc'
output_folder = '/Users/fviola/code/medslik_netcdf_utils/fusion_track_beach/output/'
coast_file = '/Users/fviola/code/medslik_coastlines/medf.txt'

# initialization
images = []


######################################################
#
# Functions needed for oil beached parcels
#
#######################################################

def dist(point,segments): # x3,y3 is the point

	x1=segments[:,0]
	x2=segments[:,2]
	y1=segments[:,1]
	y2=segments[:,3]
	x3=point[0]
	y3=point[1]
	px = x2-x1
	py = y2-y1
	something = px*px + py*py
	u =  ((x3 - x1) * px + (y3 - y1) * py) / something
	u[np.argwhere(u>1)]=1
	u[np.argwhere(u<0)]=0
	x = x1 + u * px
	y = y1 + u * py
	dx = x - x3
	dy = y - y3
	dist = np.sqrt(dx*dx + dy*dy)

	return dist

def haversine(lon1, lat1, lon2, lat2):
	"""
	Calculate the great circle distance between two points
	on the earth (specified in decimal degrees). Taken from stackoverflow, thanks to Michael Dunn!
	http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
	"""
	# convert decimal degrees to radians
	lon1=np.radians(lon1)
	lon2=np.radians(lon2)
	lat1=np.radians(lat1)
	lat2=np.radians(lat2)

	# haversine formula
	dlon = lon2 - lon1
	dlat = lat2 - lat1
	a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
	c = 2 * np.arcsin(np.sqrt(a))

	# 6367 km is the radius of the Earth
	km = 6367 * c
	return km


######################################################
#
# Function to track oil (-> 'is' = 1 or 2)
#
#######################################################

def oil_track(ncfile, ds, time_index, lats, lons):

    # load start position of the spill
    y0 = ncfile.variables['non_evaporative_volume'].initial_position_y
    x0 = ncfile.variables['non_evaporative_volume'].initial_position_x
    print ("Spill location = " + str(x0) + "W ::::: " + str(y0) + "N ")

    # variable extraction
    lats = ncfile.variables['latitude'][time_index,:]
    lons = ncfile.variables['longitude'][time_index,:]
        
    # generate output grid
    grid_min_longitude = np.min(lons)-ds
    grid_min_latitude = np.min(lats)-ds
    grid_max_longitude = np.max(lons)+ds
    grid_max_latitude = np.max(lats)+ds

    try:
        x_points = np.arange(grid_min_longitude+ds/2,grid_max_longitude,ds)
    except:
        pdb.set_trace()
    y_points = np.arange(grid_min_latitude+ds/2,grid_max_latitude,ds)
    box = np.zeros((len(y_points),len(x_points)))
    area = (ds*110)**2

    # conversion factor - barrel to tonnes
    oil_density = ncfile.variables['non_evaporative_volume'].oil_density
    #parcel_volume = ncfile.variables['non_evaporative_volume'].volume_of_parcel
    rbm3=0.158987
    barrel2tonnes=1/(rbm3*(oil_density/1000))

    # extract variables of interest
    particle_status = ncfile.variables['particle_status'][time_index,:]
    evaporative_volume = ncfile.variables['evaporative_volume'][time_index,:]
    non_evaporative_volume = ncfile.variables['non_evaporative_volume'][time_index,:]

    # oil track
    iNoise=np.logical_or(particle_status <= 0, particle_status > 2).nonzero()[0]
    lats = np.delete(lats, (iNoise), axis=0)
    lons = np.delete(lons, (iNoise), axis=0)
    evaporative_volume = np.delete(evaporative_volume, (iNoise), axis=0)
    non_evaporative_volume = np.delete(non_evaporative_volume, (iNoise), axis=0)

    xp=np.round((lons-np.min(x_points))/ds)
    yp=np.round((lats-np.min(y_points))/ds)
    total_volume=(evaporative_volume+non_evaporative_volume)/barrel2tonnes

    for aa in range(0,len(xp)):
        box[yp[aa].astype(int),xp[aa].astype(int)] = box[yp[aa].astype(int),xp[aa].astype(int)] + total_volume[aa]/area

    return x0, y0,x_points, y_points, box



######################################################
#
# Function to calculate beached oil (-> 'is' negative)
#
######################################################

def oil_beach(ncfile, ds, time_index, lats, lons):

    # variable extraction
    lats = ncfile.variables['latitude'][time_index,:]
    lons = ncfile.variables['longitude'][time_index,:]

    # extract variables of interest
    particle_status = ncfile.variables['particle_status'][time_index,:]
    iNoise=np.argwhere(particle_status >= 0)
    lats = np.delete(lats, (iNoise), axis=0)
    lons = np.delete(lons, (iNoise), axis=0)
    particle_status = np.delete(particle_status, (iNoise), axis=0)
    iBeaching=(np.transpose(lons),np.transpose(lats),np.transpose(particle_status))

    # TODO -- WIP
    # xp=np.round((lons-np.min(x_points))/ds)
    # yp=np.round((lats-np.min(y_points))/ds)
    # total_volume=(evaporative_volume+non_evaporative_volume)/barrel2tonnes

    # for aa in range(0,len(xp)):
    #     box[yp[aa].astype(int),xp[aa].astype(int)] = box[yp[aa].astype(int),xp[aa].astype(int)] + total_volume[aa]/area

    return iBeaching


######################################################
#
# main
#
######################################################

if __name__ == "__main__":

    # load MEDSLIK netCDF output file
    ncfile = Dataset(inputFile,'r')

    # prepare segment positions
    iTargetSites = np.loadtxt(coast_file)
    iSegmentLengths=haversine(iTargetSites[:,0],iTargetSites[:,1],iTargetSites[:,2],iTargetSites[:,3])
    
    # find the origin of the spill
    y0 = ncfile.variables['non_evaporative_volume'].initial_position_y
    x0 = ncfile.variables['non_evaporative_volume'].initial_position_x
    time_list = ncfile.variables['time'][:]
    print ("Spill initial location = " + str(x0) + "W ::::: " + str(y0) + "N ")

    # variable extraction
    lats = ncfile.variables['latitude'][0,:]
    lons = ncfile.variables['longitude'][0,:]
        
    # generate output grid
    grid_resolution = 0.15/110
    max_ds = np.max(time_list*.8)/110.
    grid_min_longitude = np.min(lons) - max_ds #grid_resolution
    grid_min_latitude = np.min(lats) - max_ds #grid_resolution
    grid_max_longitude = np.max(lons) + max_ds #grid_resolution
    grid_max_latitude = np.max(lats) + max_ds #grid_resolution

    # initiate coastline
    print("[oil_track_mdk2] -- Initiating Basemap")
    m = Basemap(llcrnrlon=grid_min_longitude-.1,llcrnrlat=grid_min_latitude-.1,\
                urcrnrlon=grid_max_longitude+.101,urcrnrlat=grid_max_latitude+.101,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='f',projection='merc',\
                lat_0=(grid_max_latitude + grid_min_latitude)/2,\
		lon_0=(grid_max_longitude + grid_min_longitude)/2)

    for ii in range(len(time_list)-50):

        print("[oil_track_mdk2] -- Timestep: %s" % ii)
        
        ###############################################
        #
        # map customization
        #
        ###############################################
        
        # Plot trajectory - set map
        plt.figure(figsize=(7, 7)) # This increases resolution
        
        # Plot coastline
        m.drawcoastlines(linewidth=0.05)
        m.fillcontinents(alpha=0.1)
        
        m.drawmeridians(np.arange(grid_min_longitude-.1,grid_max_longitude+.101, \
        	                  ((grid_max_longitude+.101)-(grid_min_longitude-.1))/3),\
                        labels=[0,0,0,1],color='black',linewidth=0.03) # draw parallels
        m.drawparallels(np.arange(grid_min_latitude-.1,grid_max_latitude+.101, \
                                  ((grid_max_latitude+.101)-(grid_min_latitude-.1))/3), \
                        labels=[1,0,0,0],color='black',linewidth=0.03) # draw meridians

        
        ###############################################
        #
        # plot oil track
        #
        ###############################################
        
        oil_x0, oil_y0, oil_x_points, oil_y_points, oil_box = oil_track(ncfile, grid_resolution, ii, lats, lons)
        oiled_grid_points = np.argwhere(oil_box!=0)
        X = oil_x_points[oiled_grid_points[:,1]]
        Y = oil_y_points[oiled_grid_points[:,0]]
        x, y = m(X,Y)
        box_plot = oil_box[oiled_grid_points[:,0],oiled_grid_points[:,1]]    

        # Plot trajectory
        m.scatter(x, y, s=[1.], c=np.log(box_plot), vmin=np.log(0.05),vmax=np.log(1.0), edgecolor='')

        # Plot initial spill point
        mx0, my0 = m(x0, y0)
        m.plot(mx0, my0,'k+',markersize=5, marker="x")


        ###############################################
        #
        # plot oil beach
        #
        ###############################################

        lats = ncfile.variables['latitude'][ii,:]
        lons = ncfile.variables['longitude'][ii,:]
        particle_status = ncfile.variables['particle_status'][ii,:]
        iNoise=np.argwhere(particle_status >= 0)
        lats = np.delete(lats, (iNoise), axis=0)
        lons = np.delete(lons, (iNoise), axis=0)
        particle_status = np.delete(particle_status, (iNoise), axis=0)
        iBeaching=(np.transpose(lons),np.transpose(lats),np.transpose(particle_status))

        # output matrices
        iAssignedSegment=np.zeros(np.shape(lons)[0])
        iConcentrationsParcels=np.zeros(len(iTargetSites))
        iCP=np.zeros(len(iTargetSites))

        oil_density = ncfile.variables['non_evaporative_volume'].oil_density
        parcel_volume = ncfile.variables['non_evaporative_volume'].volume_of_parcel
        rbm3=0.158987
        barrel2tonnes=1/(rbm3*(oil_density/1000))
        
        print ('...assigning parcels to coastal segments...') 
        # assign a target site to the beached parcels
        for jj in range(0,np.shape(lons)[0]):
            try:
                iParcelPosition=(iBeaching[0][jj],iBeaching[1][jj])
                iParcelDistance=dist(iParcelPosition,iTargetSites)
                iParcelDistance[np.isnan(iParcelDistance)]=9999 # border segments are removed
                iClosestSegmentDist=np.min(iParcelDistance)
            except:
                print(traceback.print_exc())
                pdb.set_trace()
                
            if iClosestSegmentDist<.1/110:
                if len(np.argwhere(iParcelDistance==iClosestSegmentDist))>1:
                    iAssignedSegment[jj]=np.argwhere(iParcelDistance==iClosestSegmentDist)[0]
                    
                else:
                    iAssignedSegment[jj]=np.argwhere(iParcelDistance==iClosestSegmentDist)
                    
            iObservedOil=((-iBeaching[2][jj])/iSegmentLengths[int(iAssignedSegment[jj])])/barrel2tonnes
            iCP[int(iAssignedSegment[jj])]=iObservedOil

        
        if len(iCP>0):
            aa=np.argwhere(iCP>0)
            x_m=(iTargetSites[:,0]+iTargetSites[:,2])/2
            y_m=(iTargetSites[:,1]+iTargetSites[:,3])/2

            x_m,y_m=m(x_m,y_m)
            sorted_concs=np.argsort(iCP)
   
            for ss in sorted_concs:
                ax1=plt.subplot(111)
            
                if iCP[ss]>0:
                    cs = plt.scatter(x_m[ss],y_m[ss], s=(iCP[ss]/np.max(iCP))*80, c = iCP[ss], vmin = np.percentile(iCP[aa],5), vmax = np.percentile(iCP[aa],95),edgecolor='',alpha=0.3,cmap='gist_rainbow')                
                    cbar = m.colorbar(cs,location='bottom',pad=0.2)#pad="5%")
                    cbar.set_label('tons/km')
                
        ###############################################
        #
        # determine oiled locations
        #
        ###############################################
        
        oil_x0, oil_y0, oil_x_points, oil_y_points, oil_box = oil_track(ncfile, grid_resolution, ii, lats, lons)
        oiled_grid_points = np.argwhere(oil_box!=0)
        X = oil_x_points[oiled_grid_points[:,1]]
        Y = oil_y_points[oiled_grid_points[:,0]]
        x, y = m(X,Y)
        box_plot = oil_box[oiled_grid_points[:,0],oiled_grid_points[:,1]]    


        ###############################################
        #
        # map customization
        #
        ###############################################
        
        #         # Colorbar setup
        #         ticks = numpy.log([0.05, 0.1, 0.5, 1.0])
        #         cbar = plt.colorbar(ticks=ticks, format='$%.2f$', orientation='horizontal')
        #         cbar.ax.set_xticklabels(['0.01','0.05','0.1','0.5','1.0'])        
        #         cbar.set_label('tons/km2')

        ###############################################
        #
        # Save and close the image
        #
        ###############################################

        plt.title('Surface oil concentrations for '+ '%03d' % (ii+1) + 'h')
        plt.savefig(output_folder + '/surface_oil_' + '%03d' % (ii+1) + 'h.png',dpi=600,bbox_inches='tight')
        plt.close('all')
        
        # add image to list
        images.append(imageio.imread(output_folder + '/surface_oil_' + '%03d' % (ii+1) + 'h.png'))

# generate output gif
print("Generating output gif in %s/anim.gif" % output_folder)
imageio.mimsave(output_folder + '/anim.gif', images, fps=2)
