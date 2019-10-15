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
	filefolder ='/location_of_directory/'

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

	np.save(opj(filefolder, 'Planet_'+str(planet_id)+'_timelines'),timelines)


def call_horizons(verbose=False, time_interval='1h'):

    Epochs={'start':"2018-10-01", 'stop':"2018-10-10", 'step':time_interval}
    obj=Horizons(id='5', id_type='majorbody', location=None, epochs=Epochs)
    eph=obj.ephemerides()

    # Converting relevant information into numpy arrays
    ra = np.array(eph['RA'])
    dec = np.array(eph['DEC'])
    t = np.array(eph['datetime_jd'])

    # Combining information into a dictionary
    data = {'ra':ra, 'dec':dec, 't':t}

    # Writing data to a pickle file
    with open('horizons_output.pkl', 'wb') as handle:
        pickle.dump(data, handle)

    # Printing out some useful debutting statements
    if verbose:

        print(eph['datetime_str'], eph['RA'], eph['DEC'])        
        print(type(eph['RA'][0]))
        print(type(eph['RA']))
        print(eph.columns)
        print(np.shape(ra))
        print(np.shape(dec))
        print(np.shape(t))

    return

def main():

    call_horizons(verbose=True)

if __name__ == '__main__':

    main()
