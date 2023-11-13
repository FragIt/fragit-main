"""
Copyright (C) 2011-2023 Casper Steinmann
"""

from fragit.template import Template


class PymolTemplate(Template):
    def __init__(self, directories, infile: str, outfile: str):
        Template.__init__(self, directories, infile, outfile)
        self._set_template_type("pymol")
        self._set_load_structure_string("cmd.do(\"load %s\")")

    def format_single_fragment(self, fragment):
        s = ""
        for atom in fragment:
            s += "%s," % (atom)
        return s[:-1] + ":"

    def format_fragments(self):
        fragments = "fragments_data=\"%s\"\n"
        s = ""
        for fragment in self.fragments_data:
            s += self.format_single_fragment(fragment)
        fragments = fragments % (s[:-1])
        return fragments

    def format_buffer(self):
        fragments = "buffer_data=\"%s\"\n"
        s = ""
        for i,fragment in enumerate(self.fragments_data):
            if self.buffer_data[i] == 2:
                s += self.format_single_fragment(fragment)
        fragments = fragments % (s[:-1])
        return fragments

    def format_active(self):
        fragments = "active_data=\"%s\"\n"
        s = fragments % (self.format_single_fragment(self.active_data)[:-1])
        return s

    def format_backbone(self):
        fragments = "backbone_data=\"%s\""
        s = fragments % (self.format_single_fragment(self.backbone_data)[:-1])
        return s

    def format_break_points(self):
        return ""

    def format_fragment_charges(self):
        charges = "fragment_charges=\"%s\""
        s = ",".join(map(str, self.fragment_charges))
        return charges % s
