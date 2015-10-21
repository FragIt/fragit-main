# FragIt

[![GitHub release](https://img.shields.io/github/release/FragIt/fragit-main.svg?style=flat)](https://github.com/FragIt/fragit-main/releases)
[![Build Status](https://travis-ci.org/FragIt/fragit-main.svg?branch=master)](https://travis-ci.org/FragIt/fragit-main)

FragIt is a python based tool that allows you to quickly fragment ["any"](http://openbabel.org/docs/2.3.0/FileFormats/Overview.html) molecule and use the produced output file(s) as an a starting point for input files in quantum chemistry programs that supports such fragment based methods.

FragIt was made out of the need to quickly benchmark many different molecules while developing new fragment based methods and is now being released in the hope that it is useful for others in their research. You can read about FragIt in the [published paper](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0044480).

Currently, FragIt supports the [fragment molecular orbital](http://en.wikipedia.org/wiki/Fragment_Molecular_Orbital) method in [GAMESS](http://www.msg.ameslab.gov/gamess/index.html). FragIt also comes with a standard purpose XYZ writer that dumps each fragment in a separate `.xyz` file. Finally, there is a new XYZ-MFCC writer to support molecular fragmentation with conjugate caps that will yield capped fragments and caps for MFCC. New output writers can be added easily to support other methods and programs.

There is also an [online-version available](http://www.fragit.org/) which requires only a browser and Java but is limited in the number of fragments / residues that can be fragmented.

## Obtaining FragIt

Since you found this file, it is obvious that you also found the source code. You can obtain the latest version from [github](https://www.github.com/FragIt/fragit-main) where [tagged releases](https://github.com/FragIt/fragit-main/releases) are also available.

## Installing FragIt

FragIt is a python library and installation is quite straight forward

    python setup.py build
    python setup.py install

to install it in the default locations. To install it in a custom location, you can run the following

    python setup.py build
    python setup.py install --prefix=/path/to/custom/installation

## Running FragIt

in a terminal, you can type

    fragit

to see its help message.

See the [wiki](https://www.github.com/FragIt/fragit-main/wiki) for more examples of how to use FragIt.

## Requirements

In order to run FragIt, you *need* the following installed on your system:

* The FragIt source code, look above for information on how to obtain it
* [Open Babel](http://www.openbabel.org) 2.3 or newer with [python language bindings](http://openbabel.org/docs/dev/Installation/install.html#compile-language-bindings) enabled. Github user [andersx](https://github.com/andersx) [wrote a guide](http://combichem.blogspot.dk/2013/12/compiling-open-babel-with-python.html) to how that is accomplished.
* [Numpy](http://numpy.scipy.org) 1.5 or newer.
* [Python](http://www.python.org) version 2.4 or later (not 3.X)
