from dataclasses import dataclass
from typing import Optional

import astropy.units as u
from astropy.coordinates import EarthLocation
from tendrils import api


@dataclass
class Site:
    siteid: int
    sitename: str
    longitude: float
    latitude: float
    elevation: float
    earth_location: Optional[EarthLocation] = None
    site_keyword: Optional[str] = None

    def __post_init__(self):
        if self.earth_location is not None:
            return
        self.earth_location = EarthLocation(
            lat=self.latitude*u.deg,
            lon=self.longitude*u.deg,
            height=self.elevation*u.m
        )

    @classmethod
    def from_flows(cls, siteid: int) -> 'Site':
        site_dict = api.get_site(siteid)
        site_dict['earth_location'] = site_dict.pop('EarthLocation')
        return cls(**site_dict)

    @classmethod
    def from_astropy(cls, sitename: str) -> 'Site':
        loc = EarthLocation.of_site(sitename)
        return cls(siteid=999, sitename=sitename,
            longitude=loc.lon.value, latitude=loc.lat.value,
            elevation=loc.height.value, earth_location=loc)

    @classmethod
    def from_query(cls) -> 'Site':
        sitename = input('Enter a site name for logging: ')
        longitude = float(input('Enter longitude in degrees: '))
        lat = float(input('Enter latitude in degrees: '))
        elevation = float(input('Enter elevation in meters: '))
        siteid = 999  # hardcoded for user defined site
        return cls(siteid=siteid, sitename=sitename,
            longitude=longitude, latitude=lat,
            elevation=elevation)
