//
//  Copyright (C) 2004-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Modified version of Embedder.h for initial conformer generation
// Using installed RDKit headers instead of local development headers

#ifndef INITIAL_EMBEDDER_H
#define INITIAL_EMBEDDER_H

#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ROMol.h>
#include <RDGeneral/types.h>

namespace RDKit {
namespace InitialEmbedder {

extern const DGeomHelpers::EmbedParameters KDG;
extern const DGeomHelpers::EmbedParameters ETDG;
extern const DGeomHelpers::EmbedParameters ETDGv2;
extern const DGeomHelpers::EmbedParameters ETKDG;
extern const DGeomHelpers::EmbedParameters ETKDGv2;
extern const DGeomHelpers::EmbedParameters ETKDGv3;
extern const DGeomHelpers::EmbedParameters srETKDGv3;

// Main function for embedding multiple conformers
// This is a standalone version that uses RDKit's installed headers
void EmbedMultipleConfs(ROMol &mol, INT_VECT &res, unsigned int numConfs,
                        DGeomHelpers::EmbedParameters &params);

}  // end of namespace InitialEmbedder
}  // end of namespace RDKit

#endif  // INITIAL_EMBEDDER_H
