#!d:\learningcentre\python\flask\aircraftmodelswebservices\openaerostruct\oas_v1\venv\scripts\python.exe
# EASY-INSTALL-ENTRY-SCRIPT: 'openmdao==3.16.0','console_scripts','openmdao'
__requires__ = 'openmdao==3.16.0'
import re
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(
        load_entry_point('openmdao==3.16.0', 'console_scripts', 'openmdao')()
    )
