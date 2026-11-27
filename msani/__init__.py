from importlib.metadata import metadata as _meta, PackageNotFoundError

try:
    _m = _meta("MolSanitizer")
    __version__ = _m["Version"]
    __author__  = _m["Author-email"]
    # SPDX metadata (PEP 639) stores the license under License-Expression.
    # A missing key raises KeyError from Python 3.15 on, so use get().
    __license__ = _m.get("License-Expression") or _m.get("License")
except PackageNotFoundError:
    # Running directly from source without installing
    __version__ = "0.0.0+dev"
    __author__  = "Thua-Phong Lam, Szymon Pach, Israel Cabeza de Vaca"
    __license__ = "Apache-2.0"

from msani.moltransform.ionizer import Ionizer
from msani.moltransform.neutralizer import Neutralizer
from msani.moltransform.tautomerizer import Tautomerizer
from msani.filtering.filters import Filters
from msani.api import Msani