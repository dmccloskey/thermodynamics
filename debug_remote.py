#VSCode remote debugging
######
import ptvsd
ptvsd.enable_attach("my_secret", address = ('0.0.0.0', 3000))
#enable the below line of code only if you want the application to wait until the debugger has attached to it
ptvsd.wait_for_attach()
######

# if the thermodynamics or component-contribution package is not located in lib/site-packages
# define the location
import sys
sys.path.append('/home/user/component-contribution')
sys.path.append('/home/user/io_utilities')
sys.path.append('/home/user/cobra_utilities')

# import the example file and run
# from thermodynamics.thermodynamics_examples import aerobicAnaerobic01
# aerobicAnaerobic01._main_()
from thermodynamics.thermodynamics_examples import aerobicAnaerobic02
aerobicAnaerobic02._main_()