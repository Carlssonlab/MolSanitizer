#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ROMol.h>
#include <memory>

int main() {
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol("CCO"));
    return mol && mol->getNumAtoms() == 3 ? 0 : 1;
}
