"""
**********************************************************************
template.py

Copyright (C) 2011-2012 Casper Steinmann

This file is part of the FragIt project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
"""

import sys
import string
import os

import util

class Template(object):
	def __init__(self, infile, outfile):
		self.infile = infile
		self.outfile = outfile
		self.data = list() # holds the template
		self.fragments_data = list()
		self.buffer_data = list()
		self.active_data = list()
		self.backbone_data = list()
		self.template_type = None
		self.load_structure_string = None

	def _setTemplateType(self,value):
		if not util.is_string(value): raise ValueError("Template type is a string value.")
		self.template_type = value

	def _setLoadStructureString(self,value):
		if not util.is_string(value): raise ValueError("Load structure string is a string value.")
		self.load_structure_string = value + "\n"

	def setFragmentsData(self, value):
		self.fragments_data = value

	def setBufferData(self, value):
		self.buffer_data = value

	def setActiveData(self, value):
		self.active_data = value

	def setBackboneData(self, value):
		self.backbone_data = value

	def setPairData(self,value):
		self.fragpairs = value

	def formatSingleFragment(self, fragment_data):
		raise NotImplementedError

	def formatFragments(self, fragments):
		raise NotImplementedError

	def formatBackbone(self):
		raise NotImplementedError

	def formatBuffer(self):
		raise NotImplementedError

	def formatActive(self):
		raise NotImplementedError

	def formatBreakPoints(self):
		raise NotImplementedError

	def override(self):
		self._readTemplateFile()
		if self.load_structure_string is None: raise ValueError("Load structure string not defined")
		for i,line in enumerate(self.data):
			if "<<fragments>>" in line:
				self.data.insert(i,self.formatFragments())
				self.data.remove(line)
			if "<<buffer>>" in line:
				self.data.insert(i,self.formatBuffer())
				self.data.remove(line)
			if "<<active>>" in line:
				self.data.insert(i,self.formatActive())
				self.data.remove(line)
			if "<<backbone>>" in line:
				self.data.insert(i,self.formatBackbone())
				self.data.remove(line)

			if "<<breakpoints>>" in line:
				self.data.insert(i,self.formatBreakPoints())
				self.data.remove(line)

			if "<<pdbfilename>>" in line:
				self.data.insert(i,self.load_structure_string % (self.infile))
				self.data.remove(line)

	def _readTemplateFile(self):
		filename = self._getTemplateFilename()
		template = open(filename,'r')
		self.data = template.readlines()

	def _getTemplateFilename(self):
		if self.template_type is None: raise ValueError("Template type not set.")
		template_file = "%s_template" % (self.template_type)
		template_path = os.path.dirname(util.__file__)
		return os.path.join(template_path,template_file)

	def write(self):
		if self.template_type is None: raise ValueError("Template type not set.")
		f = open("%s.%s" % (self.outfile, self.template_type), 'w')
		f.writelines(self.data)
		f.close()
