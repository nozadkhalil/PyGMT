#!/usr/bin/env python


import os
import sys
import shutil
import numpy as np
sys.path.append('path to directory that this files are saved')
import Converts


def curvedText(txt, origin, r, ang):
	"""
	calculate rotation angle and location of the text in map.
	"""
	txtDeg = np.linspace(ang[0],ang[1],len(txt.split()))
	y = r*np.sin(np.deg2rad(txtDeg)); x = r*np.cos(np.deg2rad(txtDeg)); 
	Clon, Clat, Ch = Converts.local2llh([x, y], [origin[0], origin[1]])
	return Clon, Clat, txtDeg





