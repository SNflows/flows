from dataclasses import dataclass
from astropy.coordinates import EarthLocation
import astropy.units as u
from typing import Optional
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
        if self.earth_location is None:
            return
        self.earth_location = EarthLocation(
            lat=self.latitude*u.deg,
            lon=self.longitude*u.deg,
            height=self.elevation*u.m
        )

    @classmethod
    def from_flows(cls, siteid: int) -> 'Site':
        return cls(**api.get_site(siteid))
    
    @classmethod
    def from_astropy(cls, sitename: str) -> 'Site':
        loc = EarthLocation.of_site(sitename)
        return cls(siteid=999, sitename=sitename, 
            longitude=loc.long.value, latitude=loc.lat.value,
            elevation=loc.height.value, earth_location=loc)
        
    @classmethod
    def from_query(cls) -> 'Site':
        sitename = input('Enter a site name for logging: ')
        longitude = float(input('Enter longitude in degrees: '))
        lat = float(input('Enter latitude in degrees: '))
        elevation = float(input('Enter elevation in meters: '))
        siteid = 1000  # hardcoded for user defined site
        return cls(siteid=999, sitename=sitename, 
            longitude=longitude, latitude=lat,
            elevation=elevation)
        