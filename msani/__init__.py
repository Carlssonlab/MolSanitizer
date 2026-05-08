from importlib.metadata import metadata as _meta, PackageNotFoundError

try:
    _m = _meta("MolSanitizer")
    __version__ = _m["Version"]
    __author__  = _m["Author-email"]
    __license__ = _m["License"]
except PackageNotFoundError:
    # Running directly from source without installing
    __version__ = "0.0.0+dev"
    __author__  = "Thua-Phong Lam, Szymon Pach, Israel Cabeza de Vaca"
    __license__ = "GPLv2"

from msani.moltransform.ionizer import Ionizer
from msani.moltransform.neutralizer import Neutralizer
from msani.moltransform.tautomerizer import Tautomerizer
from msani.filtering.filters import Filters
from msani.api import Msani