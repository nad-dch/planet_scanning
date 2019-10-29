import os
import datetime
import numpy as np
import pickle
import astroquery
import astropy
import datetime
from astroquery.jplhorizons import Horizons
opj = os.path.join

def get_planet_timelines(planet_id, start_date, end_date,
						time_step):

    #planet_id can be either a 3digit number from the IAU
    #catalogue either the barycenter of the planet
    #insert date in form YYYYMMDD
    filefolder ='//Users/nadia/Downloads/'

    #convert the dates into the correct form YYYY-MM-DD
    start_date = str(start_date)[:4]+'-'+str(start_date)[4:6]+'-'+str(start_date)[6:]
    end_date = str(end_date)[:4]+'-'+str(end_date)[4:6]+'-'+str(end_date)[6:]
    time_step = str(time_step)+'d'

    Epochs = {'start':start_date,'stop':end_date,
    'step':time_step}

    obj = Horizons(id = str(planet_id), id_type = 'majorbody',
    location = None, epochs = Epochs)

    eph = obj.ephemerides()
    ra = np.array(eph['RA'])
    dec = np.array(eph['DEC'])
    t = np.array(eph['datetime_jd'])
    timelines = {'t':t, 'ra':ra, 'dec':dec}

    np.save(opj(filefolder, 'Planet_'+str(planet_id)+'_timelines'), timelines)
    return(ra,dec)


