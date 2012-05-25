"""
**********************************************************************
pymol.py

Copyright (C) 2011-2012 Casper Steinmann

This file is part of the FragIt project.

FragIt is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FragIt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
***********************************************************************/
"""

from template import Template

class PymolTemplate(Template):
	def __init__(self,infile,outfile):
		Template.__init__(self,infile,outfile)
		self._setTemplateType("pymol")
		self._setLoadStructureString("cmd.do(\"load %s\")")

	def formatSingleFragment(self, fragment):
		s = ""
		for atom in fragment:
			s+="%s," % (atom)
		return s[:-1] + ":"

	def formatFragments(self):
		fragments = "fragments_data=\"%s\"\n"
		s = ""
		for fragment in self.fragments_data:
			s += self.formatSingleFragment(fragment)
		fragments = fragments % (s[:-1])
		return fragments

	def formatBuffer(self):
		fragments = "buffer_data=\"%s\"\n"
		s = ""
		for i,fragment in enumerate(self.fragments_data):
			if self.buffer_data[i] == 2:
				s += self.formatSingleFragment(fragment)
		fragments = fragments % (s[:-1])
		return fragments

	def formatActive(self):
		fragments = "active_data=\"%s\"\n"
		s = fragments % (self.formatSingleFragment(self.active_data)[:-1])
		return s

	def formatBackbone(self):
		fragments = "backbone_data=\"%s\""
		s = fragments % (self.formatSingleFragment(self.backbone_data)[:-1])
		return s

	def formatBreakPoints(self):
		return ""

