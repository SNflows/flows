import os
import logging
from typing import Optional, Protocol
from configparser import ConfigParser

logger = logging.getLogger(__name__)


class DirectoryProtocol(Protocol):
    archive_local: str
    output_folder: str

    def set_output_dirs(self, target_name: str, fileid: int) -> None:
        ...

    def image_path(self, image_path: str) -> str:
        ...


class Directories:
    """
    Class for creating input and output directories, given a configparser instance.
    Overwrite archive_local or output_folder to manually place elsewhere.
    """
    archive_local: Optional[str] = None
    output_folder: Optional[str] = None

    def __init__(self, config: ConfigParser):
        self.config = config

    def set_output_dirs(self, target_name: str, fileid: int) -> None:
        """
        The function is meant to be called from within a context where a
        target_name and fileid are defined, so that the output_folder
        can be created appropriately.

        Parameters:
            target_name (str): Target name.
            fileid (int): The fileid of the file being processed
        """

        # Checking for None allows manual declarations to not be overwritten.
        if self.archive_local is None:
            self.archive_local = self._set_archive()
        if self.archive_local is None:
            self.output_folder = self._set_output(target_name, fileid)

        # Create output folder if necessary.
        os.makedirs(self.output_folder, exist_ok=True)
        logger.info("Placing output in '%s'", self.output_folder)

    def _set_archive(self) -> str:
        archive_local = self.config.get('photometry', 'archive_local', fallback=None)
        if archive_local is not None and not os.path.isdir(archive_local):
            raise FileNotFoundError("ARCHIVE is not available: " + archive_local)
        logger.info(f"Using data from: {archive_local}.")
        return archive_local

    def _set_output(self, target_name: str, fileid: int) -> str:
        """
        Directory for output, defaults to current
        directory if config is invalid or empty.
        """
        output_folder_root = self.config.get('photometry', 'output', fallback='.')
        output_folder = os.path.join(output_folder_root, target_name, f'{fileid:05d}')
        return output_folder

    def image_path(self, image_path: str) -> str:
        return os.path.join(self.archive_local, image_path)
