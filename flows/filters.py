from typing import Optional
from .utilities import create_logger
logger = create_logger()
FILTERS = {
    'up': 'u_mag',
    'gp': 'g_mag',
    'rp': 'r_mag',
    'ip': 'i_mag',
    'zp': 'z_mag',
    'B': 'B_mag',
    'V': 'V_mag',
    'J': 'J_mag',
    'H': 'H_mag',
    'K': 'K_mag',
}

FALLBACK_FILTER = 'gp'


def get_reference_filter(photfilter: str) -> str:
    """
    Translate photometric filter into table column.

    Parameters:
        photfilter (str): photometric filter corresponding to key of FILTERS
    """

    _ref_filter = FILTERS.get(photfilter, None)
    if _ref_filter is None:
        logger.warning(f"Could not find filter {photfilter} in catalogs. "
                       f"Using default {FALLBACK_FILTER} filter.")
        _ref_filter = FILTERS[FALLBACK_FILTER]
    return _ref_filter

def clean_value(value: str) -> str:
    """
    Clean value.
    """
    return value.replace(' ', '').replace('-', '').replace('.', '').replace('_', '').lower()

COMMON_FILTERS = {
    'B': 'B', 'V': 'V', 'R': 'R', 'g': 'gp', 'r': 'rp', 
    'i': 'ip', 'u': 'up', 'z': 'zp',
    'Ks': 'K', 'Hs': 'H', 'Js': 'J',
    'Bessel-B': 'B', 'Bessel-V': 'V', 'Bessell-V': 'V', 'SDSS-U': 'up',
    'SDSS-G': 'gp', 'SDSS-R': 'rp', 'SDSS-I': 'ip', 'SDSS-Z': 'zp',
    'PS1-u': 'up', 'PS1-g': 'gp', 'PS1-r': 'rp', 'PS1-i': 'ip', 'PS1-z': 'zp',
    'PS2-u': 'up', 'PS2-g': 'gp', 'PS2-r': 'rp', 'PS2-i': 'ip', 'PS2-z': 'zp',
    'PS-u': 'up', 'PS-g': 'gp', 'PS-r': 'rp', 'PS-i': 'ip', 'PS-z': 'zp',
    'Yc': 'Y', 'Jc': 'J', 'Hc': 'H', 'Kc': 'K',
    'Yo': 'Y', 'Jo': 'J', 'Ho': 'H', 'Ko': 'K',
    "J_Open": "J", "H_Open": "H", "K_Open": "K",
    "B_Open": "B", "V_Open": "V", "R_Open": "r", "I_Open": "i", "Y_Open": "Y",
    "g_Open": "gp", "r_Open": "rp", "i_Open": "ip", "z_Open": "zp", 'u_Open': 'up',
}


COMMON_FILTERS_LOWER = {clean_value(key): value for key, value in COMMON_FILTERS.items()}
              
                  
def match_header_to_filter(header_dict: dict[str,str]) -> str:
    """
    Extract flows filter from header.
    """
    bad_keys = ["", "NONE", "Clear"]
    filt = header_dict.get("FILTER")
    if filt is not None and filt not in bad_keys:
        filt = match_filter_to_flows(filt)
        if filt is not None:
            return filt    
    
    for key, value in header_dict.items():
        if "FILT" in key.upper():
            if value not in bad_keys:
                filt = match_filter_to_flows(value)
                if filt is not None:
                    return filt
        
   
    raise ValueError("Could not determine filter from header. Add FILTER keyword with a flows filter to header.")
                  

def match_filter_to_flows(header_filter: str) -> Optional[str]:
    """
    Match filter header value to flows filter.
    """
    if header_filter in FILTERS:
        return header_filter

    values = max(header_filter.lower().split(' '), header_filter.lower().split('.'), key=len)
    
    
    filters_keys_lower = [str(key).lower() for key in FILTERS.keys()]
    for value in values:
        if value in filters_keys_lower:
            if FILTERS.get(value) is not None:
                return value
            return value.upper()
         
    
    filters_keys_lower = [clean_value(str(key)) for key in COMMON_FILTERS.keys()]
    for value in values:
        clean = clean_value(value)
        if clean in filters_keys_lower:
            return COMMON_FILTERS_LOWER.get(clean) 

    return None
