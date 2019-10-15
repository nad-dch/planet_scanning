import numpy as np
import astropy
from astroquery.jplhorizons import Horizons
import os

def brightness(T , d):

	h = 6.2606014e-034
	#v= 
	B = 2*h*(v**2)*k*T/(c**2)*4*np.pi*(d**2)

	return B

def get_planet_timelines(planet_id, start_date, end_date,
						time_step):
	
	#planet_id can be either a 3digit number from the IAU
	#catalogue either the barycenter of the planet
	#insert date in form YYYYMMDD
	filefolder ='/home/nadia/git_repos/planet_scanning/'
	opj = os.path.join

	#convert the dates into the correct form YYYY-MM-DD
	start_date = str(start_date)[:4]+'-'+str(start_date)[4:6]+'-'+str(start_date)[6:]
	end_date = str(end_date)[:4]+'-'+str(end_date)[4:6]+'-'+str(end_date)[6:]
	time_step = str(time_step)+'d'

	Epochs = {'start':start_date,'stop':end_date,	
			'step':time_step}

	obj = Horizons(id = str(planet_id), id_type = 'majorbody', 
					location = None, epochs = Epochs)
	eph = obj.ephemerides()

	timelines = np.matrix((eph['datetime_str'],eph['RA'],
			eph['DEC'])).transpose()

	np.save(opj(filefolder, 'Planet_'+str(planet_id)+'_timelines'),timelines)























	