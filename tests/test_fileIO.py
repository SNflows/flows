import pytest
import os
from flows.fileio import DirectoriesDuringTest, Directories, IOManager
from flows.target import Target
from tendrils import utils


def delete_directories(directories):
    if directories.output_folder is None:
        return
    if os.path.exists(directories.output_folder):
        os.rmdir(directories.output_folder)

def test_import_fileio():
    assert True


def test_DirectoriesDuringTest():
    input_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input/')
    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output/')
    tdir = DirectoriesDuringTest(input_dir, output_dir)
    assert tdir.input_dir == input_dir
    assert tdir.output_dir == output_dir
    delete_directories(tdir)


def test_Directories():
    config = utils.load_config()
    directories = Directories(config)
    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output/')
    directories.set_output_dirs('test', 0, output_folder_root=output_dir)
    assert directories.output_folder == os.path.join(output_dir, 'test', '00000')
    assert directories.photometry_path == os.path.join(output_dir, 'test', '00000', 'photometry.ecsv')
    assert directories.save_as('test.txt') == os.path.join(output_dir, 'test', '00000', 'test.txt')
    delete_directories(directories)


def test_IOManager():
    target = Target(ra=0, dec=0, name='test', photfilter='rp')
    input_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input/')
    output_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'output/')
    tdir = DirectoriesDuringTest(input_dir, output_dir)
    datafile = {'test': {'fileid': 0}}
    iom = IOManager(target, tdir, datafile)
    assert iom.target.name == 'test'
    assert iom.target.photfilter == 'rp'
    assert iom.directories.image_path('test.fits') == os.path.join(input_dir, 'test.fits')
    assert iom.get_filter() == 'r_mag'
    delete_directories(tdir)
