import sys
import os
import subprocess
import re
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import cpu_count

if len(sys.argv) != 2:
	print("Something is wrong!",file=sys.stderr)
	sys.exit(-1)
with open(sys.argv[1],'r') as f:
	source = f.read().strip()
	lst = ['classes_for_tco.cpp','geom.cpp']
	
	for filename in lst:
		with open(filename,'r') as f2:
			source2 = f2.read()
			source = source.replace(r'#include "%s"' % filename, source2)
	
	print(source)
	# co = re.compile(
            # r'')
    # res = url_re.findall(text)



