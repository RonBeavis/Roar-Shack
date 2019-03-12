#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
# Project information at https://wiki.thegpm.org/wiki/Roar_shack
#

import xml.sax
import sys
import json
import hashlib
import datetime
import gzip
import re

class mzMLHandler(xml.sax.ContentHandler):
	def __init__(self):
		self.inTag = set()
		self.jsms = {}
		self.str = ''
		self.n = 0

	def getN(self):
		return self.n

	def addFile(self,of,mh,gz):
		self.ofile = of
		self.mhash = mh
		self.isGz = gz
		self.iStart = 0
		self.aas = []
		self.pep = {}

	def startElement(self, tag, attrs):
		if tag == 'group' and attrs['type'] == 'model':
			if 'rt' in attrs and len(attrs['rt']) > 0:
				try:
					self.jsms['rt'] = float(attrs['rt'])
				except:
					pass
			self.jsms['sc'] = float(attrs['id']) 
			self.jsms['pm'] = float(attrs['mh'])
			self.jsms['pz'] = int(attrs['z'])
		if tag == 'group' and attrs['label'] == 'fragment ion mass spectrum':
			self.inTag.add('fragment ion mass spectrum')
		if tag == 'domain' and attrs['id'].rfind('.1.1') == len(attrs['id'])-4:
			self.inTag.add('domain')
			self.pep['se'] = attrs['seq']
			self.pep['ex'] = float(attrs['expect'])
			self.pep['io'] = 0
			if 'y_ions' in attrs:
				self.pep['io'] += int(attrs['y_ions'])
			if 'b_ions' in attrs:
				self.pep['io'] += int(attrs['b_ions'])
			self.pep['be'] = attrs['pre']
			self.pep['de'] = float(attrs['delta'])
			self.pep['af'] = attrs['post']
			self.iStart = int(attrs['start'])
		if tag == 'aa' and 'domain' in self.inTag:
			aa = {}
			aa['po'] = int(attrs['at']) - self.iStart
			aa['mo'] = float(attrs['modified'])
			if 'pm' in attrs:
				aa['ma'] = attrs['pm']
				if 'id' in attrs:
					aa['id'] = attrs['id']
			self.aas.append(aa)
		if tag == 'note' and 'fragment ion mass spectrum' in self.inTag:
			self.inTag.add('note')
			self.jsms['ti'] = ''
		if tag == 'GAML:Xdata' and attrs['units'] == 'MASSTOCHARGERATIO':
			self.inTag.add('GAML:Xdata')
		if tag == 'GAML:Ydata' and attrs['units'] == 'UNKNOWN':
			self.inTag.add('GAML:Ydata')
		if tag == 'GAML:attribute' and attrs['type'] == 'M+H':
			self.inTag.add('GAML:attribute:M+H')
		if tag == 'GAML:attribute' and attrs['type'] == 'charge':
			self.inTag.add('GAML:attribute:charge')

	def characters(self, content):
		if 'note' in self.inTag:
			self.jsms['ti'] += content
		if 'GAML:Xdata' in self.inTag:
			self.str += content
		if 'GAML:Ydata' in self.inTag:
			self.str += content

	def endElement(self, tag):
		if tag == 'note' and 'note' in self.inTag:
			self.inTag.remove('note')
		if tag == 'domain' and 'domain' in self.inTag:
			self.inTag.remove('domain')
		if tag == 'GAML:Xdata' and 'GAML:Xdata' in self.inTag:
			self.inTag.remove('GAML:Xdata')
			str = re.sub('\s+',' ',self.str)
			str = str.strip()
			vs = str.split(' ')
			self.jsms['ms'] = []
			for v in vs:
				self.jsms['ms'].append(float(v))
			self.str = ''
		if tag == 'GAML:Ydata' and 'GAML:Ydata' in self.inTag:
			self.inTag.remove('GAML:Ydata')
			str = re.sub('\s+',' ',self.str)
			str = str.strip()
			vs = str.split(' ')
			self.jsms['is'] = []
			for v in vs:
				self.jsms['is'].append(float(v))
			self.str = ''

		if tag == 'group' and 'fragment ion mass spectrum' in self.inTag:
			self.jsms['lv'] = 2
			self.jsms['pm'] = (self.jsms['pm'] - self.pep['de'] - 1.007276466812)/self.jsms['pz']
			self.jsms['pm'] += 1.007276466812
			self.jsms['np'] = len(self.jsms['ms'])
			if len(self.aas) > 0:
				self.pep['as'] = self.aas
			self.jsms['pe'] = self.pep
			str = json.dumps(self.jsms)
			self.mhash.update(str.encode())
			str += '\n'
			if self.isGz == 0:
				self.ofile.write(str)
			else:
				self.ofile.write(str.encode())
			self.n += 1
#			if self.n % 1000 == 0:
#				print('.',end='',flush=True)
#			if self.n % 10000 == 0:
#				print(' %i' % (self.n),flush=True)
			self.jsms = {}
			self.inTag.remove('fragment ion mass spectrum')
			self.aas = []
			self.pep = {}

isGz = 0
try:
	fpath = sys.argv[1]
	opath = sys.argv[2]
	ofile = None
	if opath.find('.gz') == -1 or opath.rfind('.gz') != len(opath) - 3:
#		print('Creating jsms plain text file ... ')
		ofile = open(opath,'w', encoding='utf-8')
	else:
#		print('Creating jsms gzip file ... ')
		ofile = gzip.open(opath,'wb')
		isGz = 1
except:
	print('Error opening file specified on command line:\ne.g. ">jsms_to_bioml.py INPUT_FILE OUTPUT_FILE"')
	exit()
mhash = hashlib.sha256()
info = {'format':'jsms 1.0','source':fpath,'created':str(datetime.datetime.now())}
str = json.dumps(info) + '\n'
if isGz == 0:
	ofile.write(str)
else:
	ofile.write(str.encode())
mhash.update(str.strip().encode())
parser = xml.sax.make_parser()
parser.setContentHandler(mzMLHandler())
parser.getContentHandler().addFile(ofile,mhash,isGz)
try:
	parser.parse(open(fpath,"r"))
except SAXParseException as error:
	print(error)
	exit()
except:
	print('Non-SAX error detected: nothing done')
	exit()
n = parser.getContentHandler().getN()
print('\t\t%i entries' % (n))
info = {"validation" : "sha256", "value" : mhash.hexdigest()}
str = json.dumps(info) + '\n'
if isGz == 0:
	ofile.write(str)
else:
	ofile.write(str.encode())
