"""
**********************************************************************
jmol.py

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

from template import Template

class JmolTemplate(Template):
	def __init__(self,infile,outfile):
		Template.__init__(self,infile,outfile)
		self._setTemplateType("jmol")
		self._setLoadStructureString("load %s")

	def formatSingleFragment(self, fragment, i=0, color=""):
		s = "select "
		for atom in fragment:
			s+="atomno=%s, " % (atom)
		s = s[:-2]
		return s + "\ndefine frag%i selected\ncolor %s frag%i\n" % (i,color,i)

	def formatFragments(self):
		colors = ["green", "blue", "red","cyan", "magenta", "yellow"]*len(self.fragments_data)
		fragments = "%s\n"
		s = ""
		for i,fragment in enumerate(self.fragments_data):
			s += self.formatSingleFragment(fragment,i+1,colors[i])
		fragments = fragments % (s[:-1])
		return fragments

	def formatBuffer(self):
		fragments = "select %s"
		s = ""
		for i,fragment in enumerate(self.fragments_data):
			if self.buffer_data[i] == 2:
				for atom in fragment:
					s += "atomno=%s, " % (atom)
		fragments = fragments % (s[:-2])
		return fragments + "\ndefine buffer selected\ncolor blue buffer\n"

	def formatActive(self):
		fragments = "select %s"
		s = ""
		for atom in self.active_data:
			s += "atomno=%s, " % (atom)
		return fragments % s[:-2] + "\ndefine active selected\ncolor red active\n"

	def formatBreakPoints(self):
		s = ""
		format_yellow = "draw cbaa%i circle (atomno=%i) scale 1.0 diameter 1.15 color yellow translucent\n"
		format_orange = "draw cbda%i circle (atomno=%i) scale 1.0 diameter 1.15 color orange translucent\n"
		for i,pair in enumerate(self.fragpairs):
			s += format_yellow % (i,pair[0])
			s += format_orange % (i,pair[1])

		return s

	def formatBackbone(self):
		fragments = "select %s"
		s = ""
		for atom in self.backbone_data:
			s += "atomno=%s, " % (atom)
		return fragments % s[:-2] + "\ndefine backbone selected\ncolor red backbone\n"

