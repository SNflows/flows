#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import hashlib

#--------------------------------------------------------------------------------------------------
def get_filehash(fname):
	"""Calculate SHA1-hash of file."""
	buf = 65536
	s = hashlib.sha1()
	with open(fname, 'rb') as fid:
		while True:
			data = fid.read(buf)
			if not data:
				break
			s.update(data)

	sha1sum = s.hexdigest().lower()
	if len(sha1sum) != 40:
		raise Exception("Invalid file hash")
	return sha1sum
