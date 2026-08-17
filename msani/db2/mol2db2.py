"""Convert in-memory MOL2 and solvation objects to DB2 data."""

import io
from types import SimpleNamespace

from msani.db2 import clash, hierarchy
from msani.db2.hierarchy import TooBigError


def mol2db2(
    mol2data,
    solvdata,
    clashfile=None,
    disttol=0.001,
):
    """Build DB2 text directly from in-memory topology and solvation data."""
    options = SimpleNamespace(
        atomtypefile=None,
        colortablefile=None,
        clashfile=clashfile,
        tolerance=disttol,
        covalent=False,
        verbose=False,
        timeit=False,
        limitset=9999999999,
        limitconf=9999999999,
        limitcoord=9999999999,
        maxrecursiondepth=1,
    )

    while len(mol2data.inputEnergy) < mol2data.xyzCount:
        mol2data.inputEnergy.append(9999.99)
        mol2data.inputTotalStrain.append(0.0)
        mol2data.inputMaxStrain.append(0.0)
        mol2data.inputHydrogens.append(0)

    mol2data.convertDockTypes(options.atomtypefile)
    mol2data.addColors(options.colortablefile)
    clash_decider = clash.Clash(options.clashfile)

    def hierarchy_data_generator(this_mol2data, depth=1):
        try:
            yield hierarchy.Hierarchy(
                this_mol2data,
                clash_decider,
                tolerance=options.tolerance,
                verbose=options.verbose,
                timeit=options.timeit,
                limitset=options.limitset,
                limitconf=options.limitconf,
                limitcoord=options.limitcoord,
                solvdata=solvdata,
            )
        except TooBigError as limit_error:
            if depth > options.maxrecursiondepth:
                raise
            breaks = hierarchy.computeBreaks(limit_error, options)
            for part in range(breaks + 1):
                split_mol2data = this_mol2data.copy()
                first = len(this_mol2data.atomXyz) * part // (breaks + 1)
                last = len(this_mol2data.atomXyz) * (part + 1) // (breaks + 1)
                split_mol2data.keepConfsOnly(first, last)
                if split_mol2data.atomXyz:
                    yield from hierarchy_data_generator(
                        split_mol2data, depth=depth + 1,
                    )

    output = io.StringIO()
    for hierarchy_data in hierarchy_data_generator(mol2data):
        hierarchy_data.writeFile(
            fileHandle=output,
            verbose=options.verbose,
            timeit=options.timeit,
            limitset=options.limitset,
        )
    return output.getvalue()
