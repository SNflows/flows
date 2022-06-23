import os
from pathlib import Path
from typing import Optional, Protocol, Dict, Union
from configparser import ConfigParser
from bottleneck import allnan
from tendrils import api, utils
from .load_image import load_image
from .utilities import create_logger
from .target import Target
from .image import FlowsImage
from . import reference_cleaning as refclean
from .filters import get_reference_filter
logger = create_logger()


class DirectoryProtocol(Protocol):
    archive_local: str
    output_folder: str

    def set_output_dirs(self, target_name: str, fileid: int, create: bool = True) -> None:
        raise NotImplementedError

    def image_path(self, image_path: str) -> str:
        raise NotImplementedError

    @property
    def photometry_path(self) -> str:
        raise NotImplementedError

    def save_as(self, filename: str) -> str:
        raise NotImplementedError

    @property
    def log_path(self) -> str:
        raise NotImplementedError

    @classmethod
    def from_fid(cls, fid: int) -> 'DirectoryProtocol':
        raise NotImplementedError


class Directories:
    """
    Class for creating input and output directories, given a configparser instance.
    Overwrite archive_local or output_folder to manually place elsewhere.
    """
    archive_local: Optional[str] = None
    output_folder: Optional[str] = None

    def __init__(self, config: Optional[ConfigParser] = None):
        self.config = config or utils.load_config()


    def set_output_dirs(self, target_name: str, fileid: int, create: bool = True,
                        output_folder_root: Optional[str] = None) -> None:
        """
        The function is meant to be called from within a context where a
        target_name and fileid are defined, so that the output_folder
        can be created appropriately.

        Parameters:
            target_name (str): Target name.
            fileid (int): The fileid of the file being processed
            create (bool): Whether to create the output_folder if it doesn't exist.
            output_folder_root (str): Overwrite the root directory for output.
        """

        # Checking for None allows manual declarations to not be overwritten.
        if self.archive_local is None:
            self.archive_local = self._set_archive()

        self.output_folder = self._set_output(target_name, fileid, output_folder_root)

        # Create output folder if necessary.
        if create:
            os.makedirs(self.output_folder, exist_ok=True)
            logger.info("Placing output in '%s'", self.output_folder)

    def _set_archive(self) -> Optional[str]:
        archive_local = self.config.get('photometry', 'archive_local', fallback=None)
        if archive_local is not None and not os.path.isdir(archive_local):
            raise FileNotFoundError("ARCHIVE is not available: " + archive_local)
        logger.info(f"Using data from: {archive_local}.")
        return archive_local

    def _set_output(self, target_name: str, fileid: int, output_folder_root: Optional[str] = None) -> str:
        """
        Directory for output, defaults to current
        directory if config is invalid or empty.
        """
        output_folder_root = self.config.get('photometry', 'output', fallback='.') if output_folder_root is None \
            else output_folder_root
        output_folder = os.path.join(output_folder_root, target_name, f'{fileid:05d}')
        return output_folder

    def image_path(self, image_path: str) -> str:
        return os.path.join(self.archive_local, image_path)

    @property
    def photometry_path(self) -> str:
        return os.path.join(self.output_folder, 'photometry.ecsv')

    def save_as(self, filename: str) -> str:
        return os.path.join(self.output_folder, filename)

    @property
    def log_path(self) -> str:
        return os.path.join(self.output_folder, "photometry.log")

    @classmethod
    def from_fid(cls, fid: int, config: Optional[ConfigParser] = None, create: bool = True,
                 datafile: Optional[Dict] = None) -> DirectoryProtocol:
        instance = cls(config)
        datafile = datafile or api.get_datafile(fid)
        instance.set_output_dirs(datafile['target_name'], fid, create)
        return instance


class DirectoriesDuringTest:
    """
    Directory class in testing config.
    """
    archive_local = None
    output_folder = None

    def __init__(self, input_dir: Union[str, Path], output_dir: Union[str, Path]):
        self.input_dir = input_dir
        self.output_dir = output_dir

    def set_output_dirs(self, target_name: str, fileid: int, create: bool = True) -> None:
        self.output_folder = os.path.join(self.output_dir, target_name, f'{fileid:05d}')
        if create:
            os.makedirs(self.output_folder, exist_ok=True)

    def image_path(self, image_path: str) -> str:
        return os.path.join(self.input_dir+image_path)

    @property
    def photometry_path(self) -> str:
        return os.path.join(self.output_folder, 'photometry.ecsv')

    def save_as(self, filename: str) -> str:
        return os.path.join(self.output_folder, filename)

    @property
    def log_path(self) -> str:
        return os.path.join(self.output_folder, "photometry.log")

    @classmethod
    def from_fid(cls, fid: int, input_dir: Union[str, Path] = './test/input',
                 output_dir: Union[str, Path] = './test/output', datafile: Optional[Dict] = None) -> 'DirectoryProtocol':
        instance = cls(input_dir, output_dir)
        datafile = datafile or api.get_datafile(fid)
        instance.set_output_dirs(datafile['target_name'], fid)
        return instance


class IOManager:
    """
    Implement a runner to shuffle data.
    """

    def __init__(self, target: Target,
                 directories: DirectoryProtocol,
                 datafile: Dict):
        self.target = target
        self.directories = directories
        self.output_folder = directories.output_folder
        self.archive_local = directories.archive_local
        self.datafile = datafile
        self.diff_image_exists = False

    def _load_image(self, image_path: str) -> FlowsImage:
        """
        Load an image from a file.
        """
        # Load the image from the FITS file:
        image = load_image(self.directories.image_path(image_path), target_coord=self.target.coords)
        return image

    def load_science_image(self, image_path: Optional[str] = None) -> FlowsImage:
        image_path = image_path or self.datafile['path']
        image = self._load_image(image_path)
        logger.info("Load image '%s'", self.directories.image_path(image_path))
        image.fid = self.datafile['fileid']
        image.template_fid = None if self.datafile.get('template') is None else self.datafile['template']['fileid']
        return image

    def get_filter(self):
        return get_reference_filter(self.target.photfilter)

    def load_references(self, catalog: Optional[Dict] = None) -> refclean.References:
        use_filter = self.get_filter()
        references = api.get_catalog(self.target.name)['references'] if catalog is None else catalog['references']
        references.sort(use_filter)
        # Check that there actually are reference stars in that filter:
        if allnan(references[use_filter]):
            raise ValueError("No reference stars found in current photfilter.")
        return refclean.References(table=references)

    def load_diff_image(self) -> Optional[FlowsImage]:
        diffimage_df = self.datafile.get('diffimg', None)
        diffimage_path = diffimage_df.get('path', None) if diffimage_df else None
        self.diff_image_exists = diffimage_path is not None
        if diffimage_df and not self.diff_image_exists:
            logger.warning("Diff image present but without path, skipping diff image photometry")
        if self.diff_image_exists:
            diffimage = self._load_image(diffimage_path)
            logger.info("Load difference image '%s'", self.directories.image_path(diffimage_path))
            diffimage.fid = diffimage_df['fileid']
            return diffimage
        return None

    @classmethod
    def from_fid(cls, fid: int, directories: Optional[DirectoryProtocol] = None,
                 create_directories: bool = True, datafile: Optional[Dict] = None) -> 'IOManager':
        """
        Create an IOManager from a fileid.
        """
        datafile = datafile or api.get_datafile(fid)
        target = Target.from_fid(fid=fid, datafile=datafile)
        directories = directories or Directories.from_fid(fid=fid, create=create_directories, datafile=datafile)
        return cls(target=target, directories=directories, datafile=datafile)


def del_dir(target: Union[Path, str],
            only_if_empty: bool = False,
            delete_parent_if_file: bool = False) -> None:
    """
    Delete a given directory and its subdirectories.
    Ex: If filename is the path to a file: `del_dir(Path(filename).parent)`

    :param target: The directory to delete
    :param only_if_empty: Raise RuntimeError if any file is found in the tree
    :param delete_parent_if_file: Delete the parent directory if it is a file
    """
    target = Path(target).expanduser()
    if not target.exists():
        logger.warning("Not deleted: Directory '%s' does not exist", target)
        return
    if not target.is_dir() and delete_parent_if_file:
        logger.warning("Not deleted: '%s' is not a directory and delete_parent_if_file was False", target)
        return
    for p in sorted(target.glob('**/*'), reverse=True):  # This should also work on Windows (fingers crossed).
        if not p.exists():
            continue
        p.chmod(0o666)  # This should also work on Windows but we should have read/write permissions anyway
        if p.is_dir():
            p.rmdir()
        else:
            if only_if_empty:
                raise RuntimeError(f'{p.parent} is not empty!')
            p.unlink()
    target.rmdir()
