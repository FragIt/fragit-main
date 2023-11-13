"""
Copyright (C) 2011-2023 Casper Steinmann
"""

from fragit.template import Template


class JmolTemplate(Template):
    def __init__(self, directories, infile, outfile):
        Template.__init__(self, directories, infile, outfile)
        self._set_template_type("jmol")
        self._set_load_structure_string("load %s")

    def format_single_fragment(self, fragment, i=0, color=""):
        s = "select "
        for atom in fragment:
            s += "atomno=%s, " % atom
        s = s[:-2]
        return s + "\ndefine frag%i selected\ncolor %s\n" % (i, color)

    def format_fragments(self):
        colors = ["green", "blue", "red","cyan", "magenta", "yellow"]*len(self.fragments_data)
        fragments = "%s\n"
        s = ""
        for i, fragment in enumerate(self.fragments_data):
            s += self.format_single_fragment(fragment, i + 1, colors[i])
        fragments = fragments % (s[:-1])
        return fragments

    def format_buffer(self):
        fragments = "select %s"
        s = ""
        for i, fragment in enumerate(self.fragments_data):
            if self.buffer_data[i] == 2:
                for atom in fragment:
                    s += "atomno=%s, " % atom
        fragments = fragments % (s[:-2])
        return fragments + "\ndefine buffer selected\ncolor blue buffer\n"

    def format_active(self):
        fragments = "select %s"
        s = ""
        for atom in self.active_data:
            s += "atomno=%s, " % atom
        return fragments % s[:-2] + "\ndefine active selected\ncolor red active\n"

    def format_break_points(self):
        s = ""
        format_yellow = "draw cbaa%i circle (atomno=%i) scale 1.0 diameter 1.15 color yellow translucent\n"
        format_orange = "draw cbda%i circle (atomno=%i) scale 1.0 diameter 1.15 color orange translucent\n"
        for i, pair in enumerate(self.fragpairs):
            s += format_yellow % (i, pair[0])
            s += format_orange % (i, pair[1])

        return s

    def format_backbone(self):
        fragments = "select %s"
        s = ""
        for atom in self.backbone_data:
            s += "atomno=%s, " % atom
        return fragments % s[:-2] + "\ndefine backbone selected\ncolor red backbone\n"

    def format_fragment_charges(self):
        return ""
