#! /usr/bin/env python

import os
import sys
from distutils.core import setup
from Pyrex.Distutils import build_ext, Extension
from distutils import sysconfig

sysconfig.get_config_vars("CFLAGS") # Set gcc's flags
sysconfig._config_vars["CFLAGS"] = "-fno-strict-aliasing -DNDBUG -g -Wall"

bleu_path = os.path.join("..", "bleu")
khmer_path = os.path.join("..", "..", "lib")

setup(
	name = "bleu",
	ext_modules=[
		Extension(
			"bleu",
			["bleu.pyx"],
			language="c++", # [RCK] per the pyrex/cython c++ wrapper tutorial:  http://wiki.cython.org/WrappingCPlusPlus
#			extra_objects=extraObjs,
			pyrex_cplus=True, # [RCK] does this do something different than language="c++"?
			include_dirs=[bleu_path, khmer_path])
],
	cmdclass = {'build_ext': build_ext}
)
