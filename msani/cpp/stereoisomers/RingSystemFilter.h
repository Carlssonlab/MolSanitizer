//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Vendored into MolSanitizer from RDKit master (commit f2230854ee,
// rdkit/rdkit#9659, "Add ring system 3D-feasibility filter to
// EnumerateStereoisomers"), which is not part of an RDKit release yet.
// Keep the patterns in sync with upstream.
//
#pragma once

#include <set>
#include <utility>
#include <vector>

namespace RDKit {
class ROMol;
}

using FilterPattern_t =
    std::pair<RDKit::ROMol, std::vector<std::pair<unsigned int, int>>>;

const std::vector<FilterPattern_t> &getRingSystemFilterPatterns();

void getRingPatternsParityRelations(
    const RDKit::ROMol &mol,
    std::set<std::pair<unsigned int, unsigned int>> &same,
    std::set<std::pair<unsigned int, unsigned int>> &opposite);
