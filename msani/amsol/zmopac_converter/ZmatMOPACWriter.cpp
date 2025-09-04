#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <iomanip>

// Define M_PI for Windows/MSVC compatibility
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Simple 3D point structure
struct Point3D {
    double x, y, z;
    
    Point3D() : x(0), y(0), z(0) {}
    Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    Point3D operator-(const Point3D& other) const {
        return Point3D(x - other.x, y - other.y, z - other.z);
    }
    
    double dot(const Point3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    Point3D cross(const Point3D& other) const {
        return Point3D(y * other.z - z * other.y,
                      z * other.x - x * other.z,
                      x * other.y - y * other.x);
    }
    
    double magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }
    
    Point3D normalize() const {
        double mag = magnitude();
        if (mag == 0) return Point3D();
        return Point3D(x / mag, y / mag, z / mag);
    }
};

// Simple atom structure
struct Atom {
    std::string symbol;
    Point3D position;
    std::vector<int> bonds;  // indices of bonded atoms
    
    Atom() {}
    Atom(const std::string& sym, const Point3D& pos) : symbol(sym), position(pos) {}
};

// Simple molecule structure
class Molecule {
public:
    std::vector<Atom> atoms;
    std::string name;
    
    Molecule() : name("MOL") {}
    
    int getNumAtoms() const { return static_cast<int>(atoms.size()); }
    
    void addAtom(const std::string& symbol, const Point3D& position) {
        atoms.emplace_back(symbol, position);
    }
    
    void addBond(int atom1, int atom2) {
        if (atom1 >= 0 && atom1 < atoms.size() && atom2 >= 0 && atom2 < atoms.size()) {
            atoms[atom1].bonds.push_back(atom2);
            atoms[atom2].bonds.push_back(atom1);
        }
    }
    
    bool areBonded(int i, int j) const {
        if (i >= atoms.size() || j >= atoms.size()) return false;
        return std::find(atoms[i].bonds.begin(), atoms[i].bonds.end(), j) != atoms[i].bonds.end();
    }
};

// Geometry calculation functions
double calculateDistance(const Point3D& p1, const Point3D& p2) {
    Point3D diff = p1 - p2;
    return diff.magnitude();
}

double calculateAngle(const Point3D& p1, const Point3D& p2, const Point3D& p3) {
    // Calculate angle p1-p2-p3 (p2 is the vertex)
    Point3D v1 = p1 - p2;
    Point3D v2 = p3 - p2;
    
    double dot_product = v1.dot(v2);
    double mag1 = v1.magnitude();
    double mag2 = v2.magnitude();
    
    if (mag1 == 0 || mag2 == 0) return 0.0;
    
    double cos_angle = dot_product / (mag1 * mag2);
    // Clamp to [-1, 1] to avoid numerical errors
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    
    return acos(cos_angle) * 180.0 / M_PI;
}

double calculateDihedral(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4) {
    // Calculate dihedral angle p1-p2-p3-p4
    Point3D v1 = p2 - p1;
    Point3D v2 = p3 - p2;
    Point3D v3 = p4 - p3;
    
    Point3D n1 = v1.cross(v2).normalize();
    Point3D n2 = v2.cross(v3).normalize();
    
    double cos_angle = n1.dot(n2);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    
    double angle = acos(cos_angle) * 180.0 / M_PI;
    
    // Determine sign using scalar triple product
    Point3D cross_n1_n2 = n1.cross(n2);
    if (v2.dot(cross_n1_n2) < 0) {
        angle = -angle;
    }
    
    // Convert to positive angle [0, 360)
    if (angle < 0) {
        angle += 360.0;
    }
    
    return angle;
}

// File parsers
bool parseMol2File(const std::string& filename, Molecule& mol) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    bool in_atoms = false;
    bool in_bonds = false;
    int num_atoms = 0;
    int num_bonds = 0;
    
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) continue;
        
        if (line == "@<TRIPOS>MOLECULE") {
            if (std::getline(file, line)) {
                mol.name = line;
            }
            if (std::getline(file, line)) {
                std::istringstream iss(line);
                iss >> num_atoms >> num_bonds;
            }
            continue;
        }
        
        if (line == "@<TRIPOS>ATOM") {
            in_atoms = true;
            in_bonds = false;
            continue;
        }
        
        if (line == "@<TRIPOS>BOND") {
            in_atoms = false;
            in_bonds = true;
            continue;
        }
        
        if (line.substr(0, 8) == "@<TRIPOS>") {
            in_atoms = false;
            in_bonds = false;
            continue;
        }
        
        if (in_atoms) {
            std::istringstream iss(line);
            int atom_id;
            std::string atom_name, atom_type;
            double x, y, z;
            
            if (iss >> atom_id >> atom_name >> x >> y >> z >> atom_type) {
                // Extract element symbol from atom_type (e.g., "C.3" -> "C")
                std::string symbol = atom_type.substr(0, atom_type.find('.'));
                mol.addAtom(symbol, Point3D(x, y, z));
            }
        }
        
        if (in_bonds) {
            std::istringstream iss(line);
            int bond_id, atom1, atom2;
            std::string bond_type;
            
            if (iss >> bond_id >> atom1 >> atom2 >> bond_type) {
                // Convert from 1-based to 0-based indexing
                mol.addBond(atom1 - 1, atom2 - 1);
            }
        }
    }
    
    file.close();
    return mol.getNumAtoms() > 0;
}

bool parseXYZFile(const std::string& filename, Molecule& mol) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int num_atoms;
    
    // Read number of atoms
    if (!std::getline(file, line)) return false;
    std::istringstream(line) >> num_atoms;
    
    // Read comment line (molecule name)
    if (std::getline(file, line)) {
        mol.name = line.empty() ? "MOL" : line;
    }
    
    // Read atoms
    for (int i = 0; i < num_atoms; ++i) {
        if (!std::getline(file, line)) break;
        
        std::istringstream iss(line);
        std::string symbol;
        double x, y, z;
        
        if (iss >> symbol >> x >> y >> z) {
            mol.addAtom(symbol, Point3D(x, y, z));
        }
    }
    
    file.close();
    
    // For XYZ files, we need to guess bonds based on distances
    // Simple approach: bond if distance < sum of covalent radii * 1.2
    std::map<std::string, double> covalent_radii = {
        {"H", 0.31}, {"C", 0.76}, {"N", 0.71}, {"O", 0.66}, {"F", 0.57},
        {"P", 1.07}, {"S", 1.05}, {"Cl", 0.99}, {"Br", 1.20}, {"I", 1.39}
    };
    
    for (int i = 0; i < mol.getNumAtoms(); ++i) {
        for (int j = i + 1; j < mol.getNumAtoms(); ++j) {
            double dist = calculateDistance(mol.atoms[i].position, mol.atoms[j].position);
            
            double r1 = covalent_radii.count(mol.atoms[i].symbol) ? 
                       covalent_radii[mol.atoms[i].symbol] : 1.0;
            double r2 = covalent_radii.count(mol.atoms[j].symbol) ? 
                       covalent_radii[mol.atoms[j].symbol] : 1.0;
            
            if (dist < (r1 + r2) * 1.2) {
                mol.addBond(i, j);
            }
        }
    }
    
    return mol.getNumAtoms() > 0;
}

// Z-matrix MOPAC writer class
class ZmatMOPACWriter {
private:
    Molecule& mol;
    std::vector<std::tuple<int, int, int>> refs; // parent, grandparent, great_grandparent
    std::set<int> complete_atoms; // atoms with all three references
    
public:
    ZmatMOPACWriter(Molecule& molecule) : mol(molecule) {
        refs.resize(mol.getNumAtoms());
        computeReferences();
    }
    
private:
    int findBestParent(int atom_idx) {
        // First, try to find a complete atom that's bonded to current atom
        std::vector<int> bonded_complete;
        for (int idx : complete_atoms) {
            if (idx < atom_idx && mol.areBonded(atom_idx, idx)) {
                bonded_complete.push_back(idx);
            }
        }
        
        if (!bonded_complete.empty()) {
            return *std::max_element(bonded_complete.begin(), bonded_complete.end());
        }
        
        // If no complete atoms, find any bonded neighbor
        std::vector<int> bonded_neighbors;
        for (int idx : mol.atoms[atom_idx].bonds) {
            if (idx < atom_idx) {
                bonded_neighbors.push_back(idx);
            }
        }
        
        if (!bonded_neighbors.empty()) {
            return *std::max_element(bonded_neighbors.begin(), bonded_neighbors.end());
        }
        
        // Fallback to previous atom
        return atom_idx > 0 ? atom_idx - 1 : -1;
    }
    
    void updateAtomRefs(int atom_idx, int parent, int grandparent = -1, int great_grandparent = -1) {
        refs[atom_idx] = std::make_tuple(parent, grandparent, great_grandparent);
        
        // Add to complete set if it has all three references
        if (parent != -1 && grandparent != -1 && great_grandparent != -1) {
            complete_atoms.insert(atom_idx);
        }
    }
    
    void computeReferences() {
        int n_atoms = mol.getNumAtoms();
        
        if (n_atoms == 0) return;
        
        // Atom 0: origin
        updateAtomRefs(0, -1, -1, -1);
        
        if (n_atoms == 1) return;
        
        // Atom 1: reference to atom 0
        updateAtomRefs(1, 0, -1, -1);
        
        if (n_atoms == 2) return;
        
        // Atom 2: needs parent and grandparent
        int parent = findBestParent(2);
        if (parent == -1) parent = 1;
        
        int grandparent = -1;
        if (parent != -1) {
            grandparent = std::get<0>(refs[parent]);
        }
        
        if (grandparent == -1 || grandparent == parent) {
            // Find alternative grandparent
            for (int i = 0; i < 2; ++i) {
                if (i != parent) {
                    grandparent = i;
                    break;
                }
            }
        }
        
        updateAtomRefs(2, parent, grandparent, -1);
        
        if (n_atoms == 3) return;
        
        // Atoms 3+: need parent, grandparent, and great_grandparent
        for (int i = 3; i < n_atoms; ++i) {
            parent = findBestParent(i);
            if (parent == -1) parent = i - 1;
            
            int gparent = -1;
            int great_grandparent = -1;
            
            if (parent != -1) {
                gparent = std::get<0>(refs[parent]);
                great_grandparent = std::get<1>(refs[parent]);
            }
            
            // If parent doesn't have complete references, try alternative
            if (parent != -1 && (gparent == -1 || great_grandparent == -1)) {
                std::vector<int> bonded_complete;
                for (int idx : complete_atoms) {
                    if (idx < i && mol.areBonded(i, idx)) {
                        bonded_complete.push_back(idx);
                    }
                }
                
                if (!bonded_complete.empty()) {
                    parent = *std::max_element(bonded_complete.begin(), bonded_complete.end());
                    gparent = std::get<0>(refs[parent]);
                    great_grandparent = std::get<1>(refs[parent]);
                }
            }
            
            // Ensure all references are different and valid
            std::set<int> unique_refs = {i, parent, gparent, great_grandparent};
            unique_refs.erase(-1); // Remove invalid references
            
            if (unique_refs.size() < 4 && i >= 3) {
                // Fall back to sequential if we have duplicates
                parent = std::max(0, i - 1);
                gparent = std::max(0, i - 2);
                great_grandparent = std::max(0, i - 3);
            }
            
            // Final validation
            if (gparent == -1 && i >= 2) {
                gparent = std::max(0, i - 2);
            }
            if (great_grandparent == -1 && i >= 3) {
                great_grandparent = std::max(0, i - 3);
            }
            
            updateAtomRefs(i, parent, gparent, great_grandparent);
        }
    }
    
public:
    std::string toMOPIN() {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(6);
        
        ss << "PUT KEYWORDS HERE\n";
        ss << mol.name << "\n\n";
        
        int n_atoms = mol.getNumAtoms();
        
        for (int i = 0; i < n_atoms; ++i) {
            const Atom& atom = mol.atoms[i];
            
            int j = std::get<0>(refs[i]);
            int k = std::get<1>(refs[i]);
            int l = std::get<2>(refs[i]);
            
            if (i == 0) {
                ss << std::left << std::setw(2) << atom.symbol 
                   << std::right << std::setw(11) << 0.0 << "  1"
                   << std::setw(12) << 0.0 << "  1"
                   << std::setw(12) << 0.0 << "  1     0   0   0\n";
                   
            } else if (i == 1) {
                double d = calculateDistance(atom.position, mol.atoms[j].position);
                ss << std::left << std::setw(2) << atom.symbol 
                   << std::right << std::setw(11) << d << "  1"
                   << std::setw(12) << 0.0 << "  1"
                   << std::setw(12) << 0.0 << "  1"
                   << std::setw(6) << (j+1) << "   0   0\n";
                   
            } else if (i == 2) {
                double d = calculateDistance(atom.position, mol.atoms[j].position);
                double a = calculateAngle(atom.position, mol.atoms[j].position, mol.atoms[k].position);
                ss << std::left << std::setw(2) << atom.symbol 
                   << std::right << std::setw(11) << d << "  1"
                   << std::setw(12) << a << "  1"
                   << std::setw(12) << 0.0 << "  1"
                   << std::setw(6) << (j+1) << std::setw(4) << (k+1) << "   0\n";
                   
            } else {
                double d = calculateDistance(atom.position, mol.atoms[j].position);
                double a = calculateAngle(atom.position, mol.atoms[j].position, mol.atoms[k].position);
                double t = calculateDihedral(atom.position, mol.atoms[j].position, 
                                           mol.atoms[k].position, mol.atoms[l].position);
                ss << std::left << std::setw(2) << atom.symbol 
                   << std::right << std::setw(11) << d << "  1"
                   << std::setw(12) << a << "  1"
                   << std::setw(12) << t << "  1"
                   << std::setw(6) << (j+1) << std::setw(4) << (k+1) << std::setw(4) << (l+1) << "\n";
            }
        }
        
        return ss.str();
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_file output.mopin\n";
        std::cerr << "Supported input formats: .mol2, .xyz\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    Molecule mol;
    bool success = false;
    
    // Determine file format and parse
    if (input_file.substr(input_file.find_last_of(".") + 1) == "mol2") {
        success = parseMol2File(input_file, mol);
    } else if (input_file.substr(input_file.find_last_of(".") + 1) == "xyz") {
        success = parseXYZFile(input_file, mol);
    } else {
        std::cerr << "Error: Unsupported file format. Only .mol2 and .xyz are supported.\n";
        return 1;
    }
    
    if (!success) {
        std::cerr << "Error: Failed to parse input file " << input_file << "\n";
        return 1;
    }
    
    if (mol.getNumAtoms() == 0) {
        std::cerr << "Error: No atoms found in input file\n";
        return 1;
    }
    
    try {
        // Create writer and generate output
        ZmatMOPACWriter writer(mol);
        std::string output = writer.toMOPIN();
        
        // Write to file
        std::ofstream outfile(output_file);
        if (!outfile) {
            std::cerr << "Error: Failed to open output file " << output_file << "\n";
            return 1;
        }
        
        outfile << output;
        outfile.close();
        
        std::cout << "Successfully converted " << input_file << " to " << output_file << "\n";
        std::cout << "Processed " << mol.getNumAtoms() << " atoms\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
