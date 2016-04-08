import sys
import os
import subprocess
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import cpu_count

dirname = "seed/"
os.makedirs(dirname, exist_ok=True)

def generate_input(seed):
	filename = os.path.join(dirname,"%03d.txt" % seed)
	try:
		out_data = subprocess.check_output(
			["java","CutTheRootsVis","-seed",str(seed),"-exec",'./utils/input_saver.exe %s' %filename], timeout=10
		)
		status = "generated"
	except subprocess.TimeoutExpired:
		status = "TLE(5s)"
	except KeyboardInterrupt:
		return
	except:
		status = "RE"
	print("%03d"%seed,status)


seeds = range(1, 101)
p = Pool(processes=cpu_count())
p.map(generate_input, seeds)

	


