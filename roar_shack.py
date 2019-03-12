#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
# Project information at https://wiki.thegpm.org/wiki/Roar_shack
#

import json
import sys
import numpy
import operator
import os
import random
import time
import re

mds = 0
mod_string = {}
mod_totals = {}

def jsms_parse(_in,_mz,_del,_max,_z,_max_e):
	sp = []
	f = open(_in,'r')
	global mds
	global mod_string
	global mod_totals
	a = 0
	for l in f:
		m = re.search('"pm": (.+?)\,',l)
		if m is None:
			continue
		a += 1
		mz = float(m.group(1))
		if abs(mz - _mz) > _del:
			continue
		m = re.search('"pz": (.+?)\,',l)
		if m is None:
			continue
		z = int(m.group(1))
		if z != _z:
			continue
		s = json.loads(l)
		if 'lv' not in s:
			continue
		a += 1
		if s['pe']['ex'] > _max_e:
			continue
		seq = s['pe']['se']
		if 'as' in s['pe']:
			m = str(s['pe']['as'])+seq
			if m not in mod_string:
				mod_string[m] = mds
				mds += 1
			seq += '[%s]' % (mod_string[m])
		seq += ' ' + str(s['pz'])
		s['is'] = None
		if seq in mod_totals:
			mod_totals[seq] += 1
		else:
			mod_totals[seq] = 1
		if mod_totals[seq] < _max:
			sp.append(s)
	f.close()
	return (sp,a)

def load_unknown(_in):
	f = open(_in,'r')
	sp = []
	for l in f:
		if l.find('{') == -1:
			continue
		s = json.loads(l)
		if 'lv' not in s:
			continue
		sp.append(s)
	return sp

start = time.time()
start_full = start
ss = []
max_cluster = 100
min_cluster = 10
fp = 0
use_labels = False
z = 2
total = 10
max_expect = 0.001
resolution = 1.0
parent_mz = 0.0
parent_delta = 0.1
for v in sys.argv:
	if v.find('-h') != -1:
		print('''
>roar_shack.py lib_dir input_file (-d1.0) (-l) (-e0.001) (-M100) (-r1.0)
   where:
      lib_dir: valid directory containing .jsms spectrum libraries
   input_file: a valid .jsms file containing unknown spectra 
          -dX: parent ion mass accuracy (Da) (d = 0.1)
          -eX: maximum E-value allowed (e = 1.0e-03)
           -l:  use labels in final plot (l = no labels)
           -MX: maximum cluster size (M = 100)
           -rX: MS/MS spectrum resolution (Da) (r = 1.0)
''')
		exit()
	if v.find('-z') == 0:
		z = int(v[2:])
	if v.find('-n') == 0:
		total = int(v[2:])
	if v.find('-l') == 0:
		use_labels = True
	if v.find('-e') == 0:
		max_expect = float(v[2:])
	if v.find('-m') == 0:
		min_cluster = float(v[2:])
	if v.find('-M') == 0:
		max_cluster = float(v[2:])
	if v.find('-r') == 0:
		resolution = float(v[2:])
	if v.find('-p') == 0:
		parent_mz = float(v[2:])
	if v.find('-d') == 0:
		parent_delta = float(v[2:])

lib_dir = re.sub('/$','',sys.argv[1])
if not os.path.isdir(lib_dir):
	print('Could not find lib_dir "%s"' % (lib_dir))
	exit()

uu = load_unknown(sys.argv[2])
parent_mz = uu[0]['pm']

os.chdir(lib_dir)
ifiles = [f for f in os.listdir('.') if (f.find('.jsms') == len(f) - 5)]

print('            max_e = %.1e' % (max_expect))
print('     cluster size = %i - %i'% (min_cluster,max_cluster))
print('       resolution = %.3f Da'% (resolution))
if use_labels:
	print('           labels = yes')
else:
	print('           labels = no')
print('library directory = %s\n' % (lib_dir))
seqs = {}
lib_size = 0
print('Loading relavent data slice')
spectra_processed = 0
for f in ifiles:
	(ssa,a) = jsms_parse(f,parent_mz,parent_delta,max_cluster,uu[0]['pz'],max_expect)
	spectra_processed += a
	ss += ssa
	fp += 1
	print('.', end='',flush=True)
print(' loaded.')
print('  files processed = %i/%i' % (fp,len(ifiles)))
print('spectra processed = %i' % (spectra_processed))
Xp = []
yp = []
Zp = []
seqs = {}
labels = []
n = 0	
spec_start = 0
spec_end = 2000
spec_res = int((spec_end-spec_start)/resolution)

for s in ss:
	seq = s['pe']['se']
	if 'as' in s['pe']:
		m = str(s['pe']['as']) + seq
		if m not in mod_string:
			continue
		seq += '[%s]' % (mod_string[m])
	seq += ' ' + str(s['pz'])
	if mod_totals[seq] < min_cluster:
		continue
	if abs(s['pm'] - parent_mz) <= parent_delta:
		mzs = s['ms']
		spec = [0] * (spec_res + 1)
		scale = spec_res/(spec_end - spec_start)
		mm = 0
		if 'as' in s['pe']:
			for k in s['pe']['as']:
				mm += k['mo']
		i = 0
		for m in mzs:
			if m > spec_start and m < spec_end:
				pos = round((m-spec_start) * scale)
				if pos < spec_res:
					spec[pos] = 1.0
			i += 1
		Xp.append(spec)
		Zp.append([s['pm']])
		s_n = seq
		if s_n in seqs:
			yp.append(float(seqs[s_n]))
		else:
			seqs[s_n] = n
			n += 1
			yp.append(float(seqs[s_n]))
			labels.append(s_n)
	else:
		print(str(s['pm']) + ' not allowed')

for s in uu:
	seq = 'unknown'
	if abs(s['pm'] - parent_mz) <= parent_delta:
		mzs = s['ms']
		spec = [0] * (spec_res + 1)
		scale = spec_res/(spec_end - spec_start)
		mm = 0
		i = 0
		for m in mzs:
			if m > spec_start and m < spec_end:
				pos = round((m-spec_start) * scale)
				if pos < spec_res:
					spec[pos] = 1.0
			i += 1
		Xp.append(spec)
		s_n = seq
		if s_n in seqs:
			yp.append(float(seqs[s_n]))
		else:
			seqs[s_n] = n
			n += 1
			yp.append(float(seqs[s_n]))
			labels.append(s_n)

ss = None
X = numpy.array(Xp)
y = numpy.array(yp)
Z = numpy.array(Zp)
print('    training sets = %i' % (len(labels)-1))
print('        load time = %.1f sec' % ((time.time()-start)))
start = time.time()

from sklearn.manifold import TSNE
tsne = TSNE(n_components=2, random_state=0)

X_2d = tsne.fit_transform(X)
print('   modelling time = %.0f sec' % ((time.time()-start)))
print('      per unknown = %.3f sec' % ((time.time()-start)/len(uu)))
start = time.time()

target_ids = range(len(labels))

from matplotlib import pyplot as plt
plt.figure(figsize=(6, 5))
plt.suptitle('%s, unknown parent m/z = %.3f ± %.3f, z = %i' % (re.sub('.+\/','',lib_dir),uu[0]['pm'],parent_delta,uu[0]['pz']))
colors = ["#"+''.join([random.choice('789ABCDEF') for j in range(6)])
		for i in range(len(labels))]
colors[-1] = '#000000'
for i, c, label in zip(target_ids, colors, labels):
    plt.scatter(X_2d[y == i, 0], X_2d[y == i, 1], c=c, label=label)
if use_labels:
	plt.legend(labels[:40],loc=2, prop={'size': 7})
print('   rendering time = %.1f sec' % ((time.time()-start)))
print('       total time = %.1f sec' % ((time.time()-start_full)))
plt.show()

