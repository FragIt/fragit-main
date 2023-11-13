"""
Copyright (C) 2016-2023 Casper Steinmann

This program converts a GAMESS output file into
the appropriate HMO format needed by FMO.

It should in no way be considered an excellent
piece of software and should, if demonstrated
to anyone, be used as an example of abhorrent
coding standards.
"""
import sys

import numpy as np

filename = sys.argv[1]
basisname = sys.argv[2]

# data storage
localized_orbitals = []
orbital_type =[] 

# parse the log file one line at a time
#
# We do not attempt to parse anything outside the two lines
#
# "EDMISTON-RUEDENBERG ENERGY LOCALIZED ORBITALS"
# ...
# ...
# ...
# "DONE WITH ENERGY LOCALIZATION"
#
# We furthermore utilize the fact that for CH4 we will
# always have 5 occupied MOs.
with open(filename, 'r') as gamlog:
    parsing = False
    for line in gamlog:
        if "DONE WITH ENERGY LOCALIZATION" in line:
            parsing = False

        if "EDMISTON-RUEDENBERG ENERGY LOCALIZED ORBITALS" in line:
            parsing = True
            offset = 5

        # unless we are actually parsing stuff, do nothing
        if not parsing:
            continue

        # count down until the real data starts to become available
        offset -= 1
        if offset > 0:
            continue

        # and now we get to parse the data. Ignore if other than C.
        line_data = line.split()
        if line_data[1] != "C":
            continue

        # finally we can store the data we need
        # orbital type to find the first pz orbital
        # and of course the MO-coefficients.
        orbital_type.append(line_data[3])
        localized_orbitals.append(list(map(float, line_data[4:])))

loc_mo = np.transpose(localized_orbitals)

# find the first instance of a pz orbital in the
# localized data. This instance marks the important
# MO coefficients we need to put into the GAMESS-FMO
# input file later on.
idx_z = orbital_type.index("Z")
val_max = 0.0
idx_max = 0

# loop through the data one row at a time
nrow, ncol = np.shape(loc_mo)
for i in range(nrow):
    value = loc_mo[i][idx_z]
    if abs(value) > val_max:
        val_max = abs(value)
        idx_max = i

# once found we are now ready to print everything
# to stdout for now. The user of this script is
# responsible for piping the output to the correct
# directory
inactive_label = "0 1"
active_label = "1 0"
floatfmt = "{0:11.6f}"
n_col_width = 5

print("{} {} {}".format(basisname, ncol, nrow))

for i_row in range(nrow):
    # figure out the label for the current batch
    # of molecular orbitals. Only one is ever the
    # active one.
    label = inactive_label
    if i_row == idx_max:
        label = active_label

    # printing of the actual data
    row_data = loc_mo[i_row]
    s = label

    # print each number in succesion on the same line
    for i, data in enumerate(row_data):
        s += floatfmt.format(data)
        # UNLESS we are at the maximum number of columns of
        # data that we want (n_col_width) then we add a newline
        # UNLESS we are at the first line or at the line which
        # is supposed to end without a newline
        if (i+1)%n_col_width == 0 and i > 0 and i+1 != ncol:
            s += "\n   "

    print(s)
