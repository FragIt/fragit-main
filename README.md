# FragIt

[![GitHub release](https://img.shields.io/github/release/FragIt/fragit-main.svg?style=flat)](https://github.com/FragIt/fragit-main/releases)

FragIt is a python based tool that allows you to quickly fragment ["any"](http://openbabel.org/docs/2.3.0/FileFormats/Overview.html) molecule and use the produced output file(s) as an a starting point for input files in quantum chemistry programs that supports such fragment based methods.

FragIt was made out of the need to quickly benchmark many different molecules while developing new fragment based methods and is now being released in the hope that it is useful for others in their research. You can read about FragIt in the [published paper](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0044480).

Currently, FragIt supports the [fragment molecular orbital](http://en.wikipedia.org/wiki/Fragment_Molecular_Orbital) method in [GAMESS](http://www.msg.ameslab.gov/gamess/index.html). FragIt also comes with a standard purpose XYZ writer that dumps each fragment in a separate `.xyz` file. Finally, there is a new XYZ-MFCC writer to support molecular fragmentation with conjugate caps that will yield capped fragments and caps for MFCC. New output writers can be added easily to support other methods and programs.

## Obtaining FragIt

Since you found this file, it is obvious that you also found the source code. You can obtain the latest version from [github](https://www.github.com/FragIt/fragit-main) where [tagged releases](https://github.com/FragIt/fragit-main/releases) are also available.
You get the source code by cloning the repository

    git clone https://github.com/FragIt/fragit-main.git

## Installing FragIt in a Conda Environment With Pip
The easiest installation option is to use conda and pip.
We have supplied an `environment.yml` file for you to use.
Simply run

    cd fragit-main
    conda env create -f environment.yml

which creates an environment called fragit.
After activating the `fragit` environment

    conda activate fragit

you can install fragit using `pip` with the following command

    pip install .

## Installing FragIt Without Conda
Without conda, the installation of FragIt becomes slightly more
tedious and of course system dependent.
If you have root access, you can install all dependencies using
the system package manager, but if you don't, then you have to
install OpenBabel manually
(Github user [andersx](https://github.com/andersx) [wrote a guide](http://combichem.blogspot.dk/2013/12/compiling-open-babel-with-python.html) to how that is accomplished.)
Finally, the recommended installation method for FragIt is to use pip again 

    pip install . --user

but it is also possible using the old school approach (deprecated and not recommended)

    python setup.py build
    python setup.py install

remember to make sure that environment variables (`PATH` and `PYTHONPATH` at least!) are
set up correctly or else FragIt will complain that it cannot find itself.

## Running FragIt

in a terminal, you can type

    fragit

to see its help message.

See the [wiki](https://www.github.com/FragIt/fragit-main/wiki) for more examples of how to use FragIt.

## Requirements

In order to run FragIt, you *need* the following installed on your system:

* The FragIt source code, look above for information on how to obtain it
* [Open Babel](http://www.openbabel.org) 3 or newer with [python language bindings](http://openbabel.org/docs/dev/Installation/install.html#compile-language-bindings) enabled.
* [Numpy](http://numpy.scipy.org) 1.5 or newer.
* [Python](http://www.python.org) 3.9 (or later).
