from typing import Optional

import pytest

from flows.instruments import INSTRUMENTS, Instrument, Site


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


if __name__ == '__main__':
    pytest.main([__file__])
