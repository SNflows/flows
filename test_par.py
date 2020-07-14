#! /usr/bin/env python3

import subprocess
import concurrent.futures
import argparse
import multiprocessing as mp
from flows import api, load_config

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run photometry pipeline in parallel')
	parser.add_argument('--fileid', type=int, default=None)
	parser.add_argument('--targetid', type=int, default=2)
	parser.add_argument('--filter', type=str, default='all', choices=['all'])
	parser.add_argument('--parallel','-p',type=int,default=1)
	args = parser.parse_args()
	
	config = load_config()
	output_folder_root = config.get('photometry', 'output', fallback='.')
	
	if args.fileid is not None:
		# Run the specified fileid:
		fileids = [args.fileid]
	else:
		# Ask the API for a list of fileids which are yet to be processed:
		fileids = api.get_datafiles(targetid=args.targetid, filt=args.filter)
		print(fileids)


	#Parallelization
	cpu_use=args.parallel
	if cpu_use> mp.cpu_count():
		print(str(cpu_use)+' is greater than cpu count, '+ 
		'using max cpu count of '+str(mp.cpu_count)+' instead.')
		cpu_use=mp.cpu_count()
	if cpu_use<1:
		cpu_use=1
		
    
	files=fileids
	args_dict = {}
	for f in files:
		cmd = 'python run_photometry.py --fileid={} '.format(f)
		args_dict[f]=(cmd, 'shell=True')
	futures_d={}
	with concurrent.futures.ThreadPoolExecutor(max_workers=cpu_use) as executor:
		for fid in args_dict:
			args = args_dict[fid] 
			futures_d[executor.submit(subprocess.run, args, shell=True)]=fid

		for future in concurrent.futures.as_completed(futures_d):
			f=futures_d[future]
			if future.exception() is not None:
				print('%r generated an exception: %s' % (f, future.exception()))