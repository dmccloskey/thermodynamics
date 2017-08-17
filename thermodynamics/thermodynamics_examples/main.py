import thermodynamics

# if the thermodynamics or component-contribution package is not located in lib/site-packages
# define the location
import sys
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/thermodynamics')
sys.path.append('C:/Users/dmccloskey-sbrg/Documents/GitHub/component-contribution')

# import the example file and run
from thermodynamics.thermodynamics_examples import aerobicAnaerobic01
aerobicAnaerobic01._main_()
from thermodynamics.thermodynamics_examples import aerobicAnaerobic02
aerobicAnaerobic02._main_()