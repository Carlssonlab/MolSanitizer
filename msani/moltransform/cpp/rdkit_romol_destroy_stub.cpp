#ifdef _MSC_VER
#include <cstdio>
#include <cstdlib>
// Ensure the stub matches the dllimport signature expected by the linker
#ifdef RDKIT_GRAPHMOL_EXPORT
#undef RDKIT_GRAPHMOL_EXPORT
#endif
#define RDKIT_GRAPHMOL_EXPORT __declspec(dllimport)
#include <GraphMol/ROMol.h>

namespace RDKit {
void ROMol::destroy() {
  std::fputs("ERROR: RDKit::ROMol::destroy() stub called. RDKit linkage issue.\n", stderr);
  std::abort();
}
}
#endif
