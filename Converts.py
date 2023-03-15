#!/usr/bin/env python

import numpy as np
from geopy import distance


def local2llh(xy, origin):
	""" conver local coordinates to geographic coordinates
	    input: xy = [x, y, z]  
	           origin = [lon, lat, elevation]
	    output
	"""
	a = 6378137.0
	e = 0.08209443794970
	xy = np.double(xy) * 1000
	origin = np.double(origin) * np.pi / 180
	M0 = a * ((1 - e ** 2 / 4 - 3 * e ** 4 / 64 - 5 * e ** 6 / 256) * origin[1] -
		(3 * e ** 2 / 8 + 3 * e ** 4 / 32 + 45 * e ** 6 / 1024) * np.sin(2 * origin[1]) +
		(15 * e ** 4 / 256 + 45 * e ** 6 / 1024) * np.sin(4 * origin[1]) -
		(35 * e ** 6 / 3072) * np.sin(6 * origin[1]))
	z = xy[1, :] != -M0
	A = (M0 + xy[1, z]) / a
	B = xy[0, z] ** 2 / a ** 2 + A ** 2
	
	llh = np.zeros((3, xy.shape[1]))
	llh[1, z] = A
	delta = np.inf
	c = 0
	
	while np.max(np.abs(delta)) > 1e-8:
		C = np.sqrt((1 - e**2 * np.sin(llh[1, z])**2)) * np.tan(llh[1, z])
		M = a * ((1 - e**2/4 - 3*e**4/64 - 5*e**6/256) * llh[1, z] -
			(3*e**2/8 + 3*e**4/32 + 45*e**6/1024) * np.sin(2*llh[1, z]) +
			(15*e**4/256 + 45*e**6/1024) * np.sin(4*llh[1, z]) -
			(35*e**6/3072) * np.sin(6*llh[1, z]))
		Mn = 1 - e**2/4 - 3*e**4/64 - 5*e**6/256 - \
			2*(3*e**2/8 + 3*e**4/32 + 45*e**6/1024) * np.cos(2*llh[1, z]) + \
			4*(15*e**4/256 + 45*e**6/1024) * np.cos(4*llh[1, z]) - \
			6*(35*e**6/3072) * np.cos(6*llh[1, z])
		Ma = M / a
		delta = -(A * (C*Ma + 1) - Ma - 0.5*(Ma**2 + B) * C) / \
			(e**2 * np.sin(2*llh[1, z]) * (Ma**2 + B - 2*A*Ma) / (4*C) +
			(A - Ma) * (C*Mn - 2/np.sin(2*llh[1, z])) - Mn)
		llh[1, z] = llh[1, z] + delta
		c = c + 1
		if c > 100:
			raise ValueError('Convergence failure.')
	
	llh[0, z]  = (np.arcsin(xy[0, z]*np.cos(origin[0])/a) / np.sin(llh[1, z])) + origin[0]
	llh[0, ~z] = xy[0, ~z]/a + origin[0]
	llh[1, ~z] = 0
	llh = np.degrees(llh)  # convert to degrees
	return llh



def km2deg(dist_km, lon, lat):
	"""
	Convert distance in Km to spherical distance in degree.
	"""
	a=(lat, lon)
	b=(lat, lon+1)
	spherical_distance_degrees = dist_km/distance.distance(a,b).km
	return spherical_distance_degrees


def dist(lon1, lat1, lon2, lat2):
	"""
	Calculate distance betweet two points. 
	"""
	a=(lat1, lon1)
	b=(lat2, lon2)	
	distance = distance.distance(a,b).km
	return distance
