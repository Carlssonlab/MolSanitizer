#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <EnumerateStereoisomers.h>
#include <memory>
#include <stdexcept>
#include <vector>
#include <string>

namespace py = pybind11;
using namespace RDKit;
using namespace EnumerateStereoisomers;

// SMILES-based Stereoisomer Enumerator class
class SmilesStereoEnumerator {
 public:
  SmilesStereoEnumerator() = delete;
  SmilesStereoEnumerator(const std::string &smiles, 
                         const StereoEnumerationOptions &options = StereoEnumerationOptions(),
                         bool verbose = false) {
    // Parse SMILES to RDKit molecule
    std::unique_ptr<ROMol> mol(SmilesToMol(smiles));
    if (!mol) {
      throw std::runtime_error("Failed to parse SMILES: " + smiles);
    }
    
    // Create the stereoisomer enumerator
    dp_enumerator.reset(new StereoisomerEnumerator(*mol, options, verbose));
  }
  
  SmilesStereoEnumerator(const SmilesStereoEnumerator &other) = delete;
  SmilesStereoEnumerator(SmilesStereoEnumerator &&other) = delete;
  SmilesStereoEnumerator &operator=(const SmilesStereoEnumerator &other) = delete;
  SmilesStereoEnumerator &operator=(SmilesStereoEnumerator &&other) = delete;
  ~SmilesStereoEnumerator() = default;

  // Return next stereoisomer as SMILES string, or None if done (for Python)
  py::object next() {
    auto iso = dp_enumerator->next();
    if (!iso) {
      return py::none();
    }
    return py::cast(MolToSmiles(*iso));
  }

  // Internal method that returns string or empty string if done (for C++)
  std::string next_string() {
    auto iso = dp_enumerator->next();
    if (!iso) {
      return "";
    }
    return MolToSmiles(*iso);
  }

  unsigned int GetStereoisomerCount() {
    return dp_enumerator->getStereoisomerCount();
  }

 private:
  std::unique_ptr<StereoisomerEnumerator> dp_enumerator;
};

// Helper function to create StereoEnumerationOptions from a Python dictionary
StereoEnumerationOptions options_from_dict(const py::dict &params) {
  StereoEnumerationOptions options;
  
  if (params.contains("maxIsomers")) {
    options.maxIsomers = params["maxIsomers"].cast<unsigned int>();
  }
  if (params.contains("onlyUnassigned")) {
    options.onlyUnassigned = params["onlyUnassigned"].cast<bool>();
  }
  if (params.contains("onlyStereoGroups")) {
    options.onlyStereoGroups = params["onlyStereoGroups"].cast<bool>();
  }
  if (params.contains("unique")) {
    options.unique = params["unique"].cast<bool>();
  }
  if (params.contains("tryEmbedding")) {
    options.tryEmbedding = params["tryEmbedding"].cast<bool>();
  }
  if (params.contains("randomSeed")) {
    options.randomSeed = params["randomSeed"].cast<int>();
  }
  if (params.contains("timeout")) {
    options.timeout = params["timeout"].cast<double>();
  }
  
  return options;
}

// Function to enumerate all stereoisomers using a Python dictionary
std::vector<std::string> enumerate_stereoisomers(const std::string &smiles,
                                                 const py::dict &params,
                                                 bool verbose = false) {
  StereoEnumerationOptions options = options_from_dict(params);
  
  std::vector<std::string> results;
  SmilesStereoEnumerator enumerator(smiles, options, verbose);
  
  std::string next_smiles;
  while (!(next_smiles = enumerator.next_string()).empty()) {
    results.push_back(next_smiles);
  }
  
  return results;
}


// PyBind11 module definition
PYBIND11_MODULE(msani_stereoisomers, m) {
    m.doc() = "The MolSanitizer C++ Stereoisomer Enumerator\n\n"
              "This module provides high-performance stereoisomer enumeration based on the rdEnumerateStereoisomers of RDKit version >= 2025.3.4.";
    
    // Version information
    m.attr("__author__") = "Phong Lam, Uppsala University (2025)";

    // Expose StereoEnumerationOptions (kept for documentation purposes)
    std::string docString = "EnumerateStereoisomers options.\n\n"
                            "NOTE: This class is kept for reference but is not directly used.\n"
                            "Use Python dictionaries with enumerate_stereoisomers() instead.";
    py::class_<StereoEnumerationOptions>(m, "StereoEnumerationOptions", docString.c_str())
        .def(py::init<>())
        .def_readwrite("tryEmbedding", &StereoEnumerationOptions::tryEmbedding,
                       "If true, the process attempts to generate a standard RDKit distance geometry"
                       " conformation for the stereoisomer. If this fails, we assume that the stereoisomer is"
                       " non-physical and don't return it. NOTE that this is computationally expensive and is"
                       " just a heuristic that could result in stereoisomers being lost. Default=False")
        .def_readwrite("onlyUnassigned", &StereoEnumerationOptions::onlyUnassigned,
                       "If true, stereocenters which have a specified stereochemistry will not be"
                       " perturbed unless they are part of a relative stereo group. Default=True.")
        .def_readwrite("onlyStereoGroups", &StereoEnumerationOptions::onlyStereoGroups,
                       "If true, only find stereoisomers that differ at the StereoGroups associated with"
                       " the molecule. Default=False.")
        .def_readwrite("unique", &StereoEnumerationOptions::unique,
                       "If true, only stereoisomers that differ in canonical SMILES will be"
                       " returned. Default=True.")
        .def_readwrite("maxIsomers", &StereoEnumerationOptions::maxIsomers,
                       "The maximum number of isomers to yield. If the number of possible isomers"
                       " is greater than maxIsomers, a random subset will be yielded. If 0, there"
                       " is no maximum. Since every additional stereocenter doubles the number of"
                       " results (and execution time) it's important to keep an eye on this.")
        .def_readwrite("randomSeed", &StereoEnumerationOptions::randomSeed,
                       "Seed for random number generator. Default=-1 means no seed.")
        .def_readwrite("timeout", &StereoEnumerationOptions::timeout,
                       "Wall-clock timeout in seconds. 0 = no limit. Enumeration stops early"
                       " and returns whatever isomers were collected within the time window."
                       " Default=0.");

    // Expose dictionary-based enumeration function
    m.def("enumerate_stereoisomers", &enumerate_stereoisomers,
          py::arg("smiles"), 
          py::arg("params"), 
          py::arg("verbose") = false,
          "Enumerate all stereoisomers of a molecule given its SMILES and a dictionary of options.\n\n"
          "Parameters\n"
          "----------\n"
          "smiles : str\n"
          "    The SMILES string of the molecule.\n"
          "params : dict\n"
          "    Dictionary with the following optional keys:\n"
          "    - maxIsomers (int): Maximum number of isomers to yield. Default=0 (no limit)\n"
          "    - onlyUnassigned (bool): Only enumerate unassigned stereocenters. Default=True\n"
          "    - onlyStereoGroups (bool): Only enumerate StereoGroups. Default=False\n"
          "    - unique (bool): Return only unique stereoisomers. Default=True\n"
          "    - tryEmbedding (bool): Validate stereoisomers via embedding. Default=False\n"
          "    - randomSeed (int): Random seed for subset selection. Default=-1\n"
          "verbose : bool, optional\n"
          "    Enable verbose output. Default=False\n\n"
          "Returns\n"
          "-------\n"
          "list of str\n"
          "    List of SMILES strings representing the stereoisomers.");
}