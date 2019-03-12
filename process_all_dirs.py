#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
# Project information at https://wiki.thegpm.org/wiki/Roar_shack
#

import sys
import os
import re
import subprocess
import time

start = time.time()
base_dir = '/mnt/ssd1/jsms/library/'
script_dir = '/mnt/ssd1/jsms/library/scripts/'
in_dir = base_dir + 'gpm/'
#print(in_dir)
#dirs = [a for a in os.listdir(in_dir) if os.path.isdir(in_dir + a)]
#for a in dirs:
for d,b,c in os.walk(in_dir):
	if len(c) == 0:
		continue
	od = re.sub('/gpm/','/lib/',d)
	os.makedirs(od, exist_ok=True)
	print('in = %s\nout = %s' % (d,od))
	files = [f for f in os.listdir(d) if (f.find('.gz') == len(f) - 3)]
	if len(files) == 0:
		print('\tin directory is empty')
	#print(files)
	os.chdir(d)
	n = 1
	s = 0
	for f in sorted(files):
		if os.path.isdir(f):
			continue
		nf = re.sub('.gz$','',f)
		ofile = od + '/%s.jsms' % (nf)
		if os.path.isfile(ofile):
			n += 1
			s += 1
			continue
		subprocess.run(['gzip','-dk',f])
		this_start = time.time()
		print('\t%i/%i. %s' % (n,len(files),nf))
		scr = script_dir + 'library_from_bioml.py'
		subprocess.run(['python3',scr,nf,ofile])
		n += 1
		os.remove(nf)
		print('\t\t%.2f min (%.2f min)' % ((time.time()-this_start)/60,(time.time()-start)/60.0))
	if n > 0:
		n -= 1
	print('files: %i, processed: %i, skipped: %i\n' % (n,n-s,s))

