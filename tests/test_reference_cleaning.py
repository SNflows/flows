import pytest
import flows.reference_cleaning as refclean
from flows.load_image import load_image
from flows.image import FlowsImage
import numpy as np
from astropy.table import Table


# @todo parametrize this with pytest to test with a list of images
FID1 = {'fid': 9378, 'fwhm': 4.14, 'path': 'tests/input/2020aatc/SN2020aatc_K_20201213_495s.fits.gz',
        'name': "SN2020aatc", 'n_force_reject_g2d': 47, 'n_clean_with_rsq_and_get_fwhm': 16,
        'n_get_clean_references': 34}


@pytest.fixture(scope='module')
def image():
    img = load_image(FID1['path'])
    img.fid = FID1['fid']
    return img


@pytest.fixture(scope='module')
def references(image):
    references = refclean.use_sep(image)
    references.table['starid'] = np.arange(len(references.table))
    return references


@pytest.fixture(scope='module')
def cleaner(image, references):
    return refclean.ReferenceCleaner(image=image, references=references)


def test_references(references):
    assert isinstance(references, refclean.References)
    assert references.table is not None
    assert references.xy is not None
    assert all(references.xy == references.table[['x', 'y']])
    assert references.coords is None
    assert references.mask is None

    with pytest.raises(AttributeError):
        _ = references.masked


def test_cleaner(cleaner):
    assert isinstance(cleaner, refclean.ReferenceCleaner)
    assert isinstance(cleaner.image, FlowsImage)
    assert isinstance(cleaner.references, refclean.References)
    assert cleaner.image.fid == 9378


@pytest.fixture(scope='module')
def maskedref(cleaner):
    return cleaner._clean_extracted_stars()


def test_force_reject_g2d(cleaner, maskedref):
    masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs = maskedref
    assert isinstance(masked_fwhms, np.ma.MaskedArray)
    assert isinstance(masked_ref_xys, np.ma.MaskedArray)
    assert isinstance(rsq_mask, np.ndarray)
    assert isinstance(masked_rsqs, np.ma.MaskedArray)
    assert len(rsq_mask) == len(masked_rsqs)
    assert len(masked_fwhms[rsq_mask].data) == FID1['n_force_reject_g2d']


def test_clean_with_rsq_and_get_fwhm(maskedref, references):
    masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs = maskedref
    fwhm, fwhm_clean_references = refclean.clean_with_rsq_and_get_fwhm(masked_fwhms, masked_rsqs, references.table)
    assert isinstance(fwhm, float)
    assert FID1['fwhm']-0.25 < fwhm < FID1['fwhm']+0.2
    assert isinstance(fwhm_clean_references, Table)
    assert len(fwhm_clean_references) == FID1['n_clean_with_rsq_and_get_fwhm']


def test_get_clean_references(references, maskedref):
    masked_fwhms, masked_ref_xys, rsq_mask, masked_rsqs = maskedref
    ordered_cleaned_references, order_index = refclean.get_clean_references(references.table, masked_rsqs)
    assert len(order_index) == FID1['n_get_clean_references']


def test_init_guess(references):
    guess = refclean.InitGuess(references)
    assert len(guess.init_guess_full) == len(references.table)
    assert len(guess.init_guess_full) == len(references.table)
    with pytest.raises(ValueError):
        len(guess.init_guess_diff)
    assert len(guess.init_guess_references) == len(references.table) - 1
