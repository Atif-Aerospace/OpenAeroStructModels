#!d:\learningcentre\python\flask\aircraftmodelswebservices\openaerostruct\oas_v1\venv\scripts\python.exe
# EASY-INSTALL-ENTRY-SCRIPT: 'openaerostruct','console_scripts','plot_wingbox'
__requires__ = 'openaerostruct'
import re
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(
        load_entry_point('openaerostruct', 'console_scripts', 'plot_wingbox')()
    )
