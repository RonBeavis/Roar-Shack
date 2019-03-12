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

def jsms_count(_in,_z,_me):
	f = open(_in,'r')
	sqs = {}
	global mds
	global mod_string
	libsize = 0
	for l in f:
		s = json.loads(l)
		if 'lv' not in s:
			continue
		libsize += 1
		if s['pe']['ex'] > _me:
			continue
		if 'pz' not in s:
			continue
		else:
			if s['pz'] != _z:
				continue
		seq = s['pe']['se']
		if 'as' in s['pe']:
			m = str(s['pe']['as'])+seq
			if m not in mod_string:
				mod_string[m] = str(mds)
				mds += 1
			seq += '[%s]' % (mod_string[m])
		seq += ' ' + str(s['pz'])
		if seq in sqs:
			sqs[seq] += 1
		else:
			sqs[seq] = 1
	f.close()
	return (sqs,libsize)

def jsms_parse(_in,_z,_allowed,_sqs,_max,_me):
	sp = []
	f = open(_in,'r')
	global mds
	global mod_string
	global mod_totals
	for l in f:
		s = json.loads(l)
		if 'lv' not in s:
			continue
		if 'pz' not in s:
			continue
		if s['pz'] != _z:
			continue
		if s['pe']['ex'] > _me:
			continue
		seq = s['pe']['se']
		if 'as' in s['pe']:
			m = str(s['pe']['as'])+seq
			seq += '[%s]' % (mod_string[m])
		seq += ' ' + str(s['pz'])
		if seq not in _allowed:
			continue
		if _sqs[seq] > _max:
			if random.randint(0,_sqs[seq]) >= _max:
				continue
		s['is'] = None
		sp.append(s)
		if seq in mod_totals:
			mod_totals[seq] += 1
		else:
			mod_totals[seq] = 1
	f.close()
	return sp

start = time.time()
ss = []
max_ss = 2.5e6
max_cluster = 50
min_cluster = 10
fp = 0
use_labels = False
z = 2
total = 10
max_expect = 0.001
resolution = 1.0
for v in sys.argv:
	if v.find('-h') != -1:
		print('''\nUsage: >tsne_generate.py lib_dir (-z2) (-n10) (-l) (-e0.001) (-m10) (-M50) (-r0.5)
       where:
              lib_dir valid directory containing .jsms spectrum libraries
              -zX parent ion charge selection (z = 2)
              -nX number of clusters to calculate (n = 10)
              -l  use labels in final plot (l = no labels
              -eX maximum expectation value allowed (e = 0.001)
              -mX minimum cluster size (m = 10)
              -MX maximum cluster size (M = 50)
              -rX mass resolution (Da) (r = 0.5)
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

lib_dir = re.sub('/$','',sys.argv[1])
if not os.path.isdir(lib_dir):
	print('Could not find lib_dir "%s"' % (lib_dir))
	exit()

os.chdir(lib_dir)
ifiles = [f for f in os.listdir('.') if (f.find('.jsms') == len(f) - 5)]

print('                n = %i' % (total))
print('                z = %i' % (z))
print('            max_e = %.1e' % (max_expect))
print('     cluster size = %i - %i'% (min_cluster,max_cluster))
print('       resolution = %.3f Da'% (resolution))
if use_labels:
	print('           labels = yes')
else:
	print('           labels = no')
print('library directory = %s\n' % (lib_dir))
print('Prescreening library')
seqs = {}
lib_size = 0
for f in ifiles:
	(sqs,l) = jsms_count(f,z,max_expect)
	lib_size += l
	for k in sqs:
		if k not in seqs:
			seqs[k] = sqs[k]
		else:
			seqs[k] += sqs[k]
	print('.', end='',flush=True)
print(' screened.')
print('Library size: %i' % (lib_size))
sorted_x = sorted(seqs.items(), key=operator.itemgetter(1))
allowed = {}
i = 0
spectra_processed = 0
for a in reversed(sorted_x):
	if a[1] < min_cluster:
		break
	allowed[a[0]] = a[1]
	spectra_processed += a[1]
	i += 1
	if i >= total:
		break
print('Loading relavent data slice')
for f in ifiles:
	if len(ss) > max_ss:
		break
	ss += jsms_parse(f,z,allowed,seqs,max_cluster,max_expect)
	fp += 1
	print('.', end='',flush=True)
print(' loaded.')
print('files processed = %i/%i' % (fp,len(ifiles)))
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
	if seq not in allowed:
		print(seq + ' not allowed')
		continue
	if s['pz'] == z:
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
		s_n = seq + ' ' + str(allowed[seq])
		if s_n in seqs:
			yp.append(float(seqs[s_n]))
		else:
			seqs[s_n] = n
			n += 1
			yp.append(float(seqs[s_n]))
			labels.append(s_n)
	else:
		print(str(z) + ' not allowed')

ss = None
X = numpy.array(Xp)
y = numpy.array(yp)
Z = numpy.array(Zp)
print(len(X),len(y))
print(len(seqs),len(labels))
print('load time = %.2f min' % ((time.time()-start)/60))
start = time.time()
############################################################
# Fit and transform with a TSNE
from sklearn.manifold import TSNE
tsne = TSNE(n_components=2, random_state=0)

############################################################
# Project the data in 2D
X_2d = tsne.fit_transform(X)
print('modelling time = %.2f min' % ((time.time()-start)/60))
start = time.time()
############################################################
# Visualize the data
target_ids = range(len(labels))

from matplotlib import pyplot as plt
plt.figure(figsize=(6, 5))
plt.suptitle('%s, z=%s, n=%i, max=%i, min=%i' % (re.sub('.+\/','',lib_dir),z,total,max_cluster,min_cluster))
colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
		for i in range(len(labels))]
for i, c, label in zip(target_ids, colors, labels):
    plt.scatter(X_2d[y == i, 0], X_2d[y == i, 1], c=c, label=label)
if use_labels:
	plt.legend(labels[:40],loc=2, prop={'size': 7})
print('rendering time = %.2f min' % ((time.time()-start)/60))
plt.show()

