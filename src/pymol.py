"""
Copyright (C) 2011-2016 Casper Steinmann
"""

from .template import Template

class PymolTemplate(Template):
    def __init__(self,directories,infile,outfile):
        Template.__init__(self,directories,infile,outfile)
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

