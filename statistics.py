import sys
import os
import subprocess
import re
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import cpu_count

if len(sys.argv) != 2:
	print(sys.argv[0], "[execfile]", file=sys.stderr)
	sys.exit(-1)

execfile = sys.argv[1]


def run_program(seed):
	NP = NR = AVGR = -1
	all = alive = ""
	score = 0.0
	
	try:
		out_data = subprocess.check_output(
			["java","CutTheRootsVis","-seed",str(seed),"-exec",execfile], timeout=10,stderr=subprocess.STDOUT
		)
		out_data = out_data.decode('utf-8')
		out_data = out_data.replace('\r\n','\n')
		# print(out_data)
		
		java_re = re.compile(r'NP = ([0-9]+) NR = ([0-9]+) AVGR = ([0-9]+)')
		res = java_re.findall(out_data)
		if res: 
			NP,NR,AVGR = res[0]
		
		java_re = re.compile(r'Length of all roots = ([.0-9]+)')
		res = java_re.findall(out_data)
		if res:
			all = res[0]
		
		java_re = re.compile(r'Length of alive roots = ([.0-9]+)')
		res = java_re.findall(out_data)
		if res:
			alive = res[0]
		java_re = re.compile(r'Score = ([.0-9]+)')
		
		res = java_re.findall(out_data)
		if res:
			score = res[0]
		status = "AC"
		print("%03d"%seed,status,file=sys.stderr)
		return [seed,NP,NR,AVGR,all,alive,score]
	except subprocess.TimeoutExpired:
		status = "TLE(10s)"
		print("%03d"%seed,status,file=sys.stderr)
		return []
	except KeyboardInterrupt:
		return []
	except:
		status = "RE"
		print("%03d"%seed,status,file=sys.stderr)
		return []
	
seeds = range(1, 101)
p = Pool(processes=cpu_count())
results = p.map(run_program, seeds)
for res in results:
	#print(",".join(map(str,res)))
	print(",".join(map(str,res[5:])))



