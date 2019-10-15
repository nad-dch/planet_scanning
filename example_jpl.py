import astroquery
import astropy
from astroquery.jplhorizons import Horizons

Epochs={'start':"2018-10-01", 'stop':"2018-10-10", 'step':"5d"}

obj=Horizons(id='5',id_type='majorbody',location=None,
             epochs=Epochs)

eph=obj.ephemerides()

print(eph['datetime_str'],eph['RA'],eph['DEC'])

print(obj.uri)