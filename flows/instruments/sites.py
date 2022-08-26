from dataclasses import dataclass
from astropy.coordinates import EarthLocation
from typing import Optional
from tendrils import api

@dataclass
class Site:
    siteid: int
    sitename: str
    longtitude: float
    latitude: float
    elevation: float
    earth_location: EarthLocation
    site_keyword: Optional[str] = None

    @classmethod
    def get_site(cls, siteid: int) -> 'Site':
        return cls(**api.get_site(siteid))