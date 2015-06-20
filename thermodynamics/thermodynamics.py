# if the comonent-contribution package is not located in lib/site-packages
# define the location
import sys
sys.path.append('C:\\Users\\dmccloskey-sbrg\\Documents\\GitHub\\component-contribution\\component-contribution')

# import the example file and run
from thermodynamics_examples import aerobicAnaerobic01
aerobicAnaerobic01._main_();