from astropy.table import Table
import astropy.units as u
from .image import FlowsImage
import logging

logger = logging.getLogger(__name__)


class ResultsTable(Table):

    def add_column_descriptions(self):
        # Descriptions of columns:
        self['used_for_epsf'].description = 'Was object used for building ePSF?'
        self['mag'].description = 'Measured magnitude'
        self['mag'].unit = u.mag
        self['mag_error'].description = 'Error on measured magnitude'
        self['mag_error'].unit = u.mag
        self['flux_aperture'].description = 'Measured flux using aperture photometry'
        self['flux_aperture'].unit = u.count / u.second
        self['flux_aperture_error'].description = 'Error on measured flux using aperture photometry'
        self['flux_aperture_error'].unit = u.count / u.second
        self['flux_psf'].description = 'Measured flux using PSF photometry'
        self['flux_psf'].unit = u.count / u.second
        self['flux_psf_error'].description = 'Error on measured flux using PSF photometry'
        self['flux_psf_error'].unit = u.count / u.second
        self['pixel_column'].description = 'Location on image pixel columns'
        self['pixel_column'].unit = u.pixel
        self['pixel_row'].description = 'Location on image pixel rows'
        self['pixel_row'].unit = u.pixel
        self['pixel_column_psf_fit'].description = 'Measured location on image pixel columns from PSF photometry'
        self['pixel_column_psf_fit'].unit = u.pixel
        self[
            'pixel_column_psf_fit_error'].description = 'Error on measured location on image pixel columns from PSF photometry'
        self['pixel_column_psf_fit_error'].unit = u.pixel
        self['pixel_row_psf_fit'].description = 'Measured location on image pixel rows from PSF photometry'
        self['pixel_row_psf_fit'].unit = u.pixel
        self[
            'pixel_row_psf_fit_error'].description = 'Error on measured location on image pixel rows from PSF photometry'
        self['pixel_row_psf_fit_error'].unit = u.pixel

    def add_metadata(tab):
        raise NotImplementedError()
        # # Meta-data:
        # tab.meta['fileid'] = fileid
        # tab.meta['target_name'] = target_name
        # tab.meta['version'] = __version__
        # tab.meta['template'] = None if datafile.get('template') is None else datafile['template']['fileid']
        # tab.meta['diffimg'] = None if datafile.get('diffimg') is None else datafile['diffimg']['fileid']
        # tab.meta['photfilter'] = photfilter
        # tab.meta['fwhm'] = fwhm * u.pixel
        # tab.meta['pixel_scale'] = pixel_scale * u.arcsec / u.pixel
        # tab.meta['seeing'] = (fwhm * pixel_scale) * u.arcsec
        # tab.meta['obstime-bmjd'] = float(image.obstime.mjd)
        # tab.meta['zp'] = zp_bs
        # tab.meta['zp_error'] = zp_error_bs
        # tab.meta['zp_diff'] = np.abs(zp_bs - zp)
        # tab.meta['zp_error_weights'] = zp_error
        # tab.meta['head_wcs'] = head_wcs  # TODO: Are these really useful?
        # tab.meta['used_wcs'] = used_wcs  # TODO: Are these really useful?

    @classmethod
    def make_results_table(cls, ref_table: Table, apphot_tbl: Table, psfphot_tbl: Table, image: FlowsImage):
        results_table = cls(ref_table)
        if len(ref_table) - len(apphot_tbl) == 1:
            results_table.add_row(0)

        if not ResultsTable.uncertainty_column_exists(psfphot_tbl):
            psfphot_tbl['flux_unc'] = psfphot_tbl['flux_fit'] * 0.02  # Assume 2% errors
            logger.warning("Flux uncertainty not found from PSF fit, assuming 2% error.")

        results_table['flux_aperture'] = apphot_tbl['flux_aperture'] / image.exptime
        results_table['flux_aperture_error'] = apphot_tbl['flux_aperture_error'] / image.exptime
        results_table['flux_psf'] = psfphot_tbl['flux_fit'] / image.exptime
        results_table['flux_psf_error'] = psfphot_tbl['flux_unc'] / image.exptime
        results_table['pixel_column_psf_fit'] = psfphot_tbl['x_fit']
        results_table['pixel_row_psf_fit'] = psfphot_tbl['y_fit']
        results_table['pixel_column_psf_fit_error'] = psfphot_tbl['x_0_unc']
        results_table['pixel_row_psf_fit_error'] = psfphot_tbl['y_0_unc']

        return results_table

    @staticmethod
    def uncertainty_column_exists(tab):
        return "flux_unc" in tab.colnames
