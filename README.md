# EOS 1 (obsolete)

Please go to: https://github.com/jianshengfeng/EOSpec_v01 for updated code.

Image analysis code (Python) for EOS 1 open-source spectrometer. 
For more information please go to:  www.erieopen.tech

Python version 2.7.15 | Anaconda 64-bit

Required libraries:
- numpy
- matplotlib
- os, sys, datetime, warnings
- PIL (needed in *ImgAna_aligncheck.py* but not in *ImgAna_minimum.py*)

Please read **ImgAna_simplified.ipynb** first. It explains how the image analysis code works.

*ImgAna_minimum.py* can be used in two ways:
- If run it as a Python script, just follow the prompt. It will first ask you to do a calibration.
- If import it as a Python module, you will have access to the **EOS1_img** class.

*ImgAna_aligncheck.py* is pretty much the same as *ImgAna_minimum.py* except that it has implemented *check_align* function. It also provides the option of automatic tilt correction (will need the PIL library to do so).
