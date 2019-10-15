import datetime
import numpy as np
import pickle
import astroquery
import astropy
import datetime
from astroquery.jplhorizons import Horizons

def call_horizons(verbose=False):

    Epochs={'start':"2018-10-01", 'stop':"2018-10-10", 'step':"1h"}

    obj=Horizons(id='5', id_type='majorbody', location=None, epochs=Epochs)
    eph=obj.ephemerides()


    ra = np.array(eph['RA'])
    dec = np.array(eph['DEC'])
    t = np.array(eph['datetime_jd'])#.to_datetime()

    data = {'ra':ra, 'dec':dec, 't':t}

    with open('horizons_output.pkl', 'wb') as handle:
        pickle.dump(data, handle)

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