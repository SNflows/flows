
import os.path
import subprocess
import shlex
import tempfile
import sys
import shutil

if __name__ == '__main__':

	photfilter = 'ip'

	reference_image = r'C:\Users\au195407\Documents\flows_archive\SN2019yvr\templates\20170420T193430_Sloan_i_cpt1m010-fl16-20170420-0097-e91.fits'
	science_image = r'C:\Users\au195407\Documents\flows_archive\SN2019yvr\elp1m008-fa05-20200107-0439-e91.fits'

	with tempfile.TemporaryDirectory() as tmpdir:
		tmpdir = 'tmpdir_12345'
		
		shutil.copy(reference_image, tmpdir)
		shutil.copy(science_image, tmpdir)

		cmd = '"{python:s}" ../ImageMatch.py -cfg "../sample.cfg" -m "{reference_image:s}" "{science_image:s}"'.format(
			python=sys.executable,
			reference_image=os.path.basename(reference_image),
			science_image=os.path.basename(science_image)
		)
		print(cmd)
		
		cmd = shlex.split(cmd)
		proc = subprocess.Popen(cmd, cwd=tmpdir)
		stdout_data, stderr_data = proc.communicate()

		print(stderr_data)
		print(stdout_data)

#python ImageMatch.py -cfg "sample.cfg" *_g_*e91.fits -m "20170420T190451_Sloan_g_cpt1m010-fl16-20170420-0091-e91.fits"
#python ImageMatch.py -cfg "sample.cfg" *_r_*e91.fits -m "20170420T192210_Sloan_r_cpt1m010-fl16-20170420-0094-e91.fits"
#python ImageMatch.py -cfg "sample.cfg" *_i_*e91.fits -m "20170420T193430_Sloan_i_cpt1m010-fl16-20170420-0097-e91.fits"
