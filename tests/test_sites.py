import tempfile
from typing import Any, Dict, List, Optional

import pytest
from matplotlib import pyplot as plt
from tendrils import api

from flows.instruments import INSTRUMENTS, Instrument, Site
from flows.target import Target
from flows.visibility import visibility


@pytest.fixture(scope='session')
def flows_sites() -> List[Dict[str, Any]]:
    return api.get_all_sites()  # type: ignore


def test_flows_sites(siteid: int = 9):

    site = Site.from_flows(siteid)

    instrument: Optional[Instrument] = None
    for instrument_name, instrument_class in INSTRUMENTS:
        if instrument_class.siteid == siteid:
            instrument = instrument_class()
    if instrument is None:
        raise ValueError(f"Expected to find site and instrument with siteid {siteid} but found None in INSTRUMENTS")

    assert instrument.telescope == 'CA 3.5m'
    assert site.siteid == siteid
    assert site.sitename == site.sitename
    assert site.site_keyword == site.site_keyword


def test_site_from_astropy_vs_flows(sitename: str = "paranal", siteid: int = 2):
    site = Site.from_astropy(sitename)
    assert site.sitename == sitename

    flows_site = Site.from_flows(siteid)
    assert int(site.latitude) == int(flows_site.latitude)
    assert int(site.longitude) == int(flows_site.longitude)


def test_user_site(monkeypatch):

    # provided inputs
    sitename = 'test'
    longitude = 11.5
    lat = 12.5
    elevation = 1200

    # creating iterator object
    answers = iter([sitename, str(longitude), str(lat), str(elevation)])

    monkeypatch.setattr('builtins.input', lambda name: next(answers))
    site = Site.from_query()

    assert site.sitename == 'test'
    assert site.siteid == 999
    assert int(site.longitude) == int(longitude)
    assert int(site.latitude) == int(lat)
    assert int(elevation) == int(elevation)

    if site.earth_location is None:
        raise ValueError(f"Expected to find site with earth location!")
    assert int(site.earth_location.lat.value) == int(lat)


# Very basic due to being moved to flows tools
def test_site_visibility(flows_sites):
    target = Target.from_tid(8)
    with pytest.deprecated_call():
        with plt.ioff():
            ax = visibility(target, siteid=2)
            assert isinstance(ax, plt.Axes)

    with tempfile.TemporaryDirectory() as tempdir:
        with pytest.deprecated_call():
            with plt.ioff():
                plotpaths = visibility(
                    target, date="2023-01-01", output=tempdir
                )

    assert not isinstance(plotpaths, plt.Axes)
    assert len(flows_sites) == len(plotpaths)


if __name__ == '__main__':
    pytest.main([__file__])
