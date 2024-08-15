import os
import subprocess
import sys
from pathlib import Path
from openbabel import openbabel as ob
from MolSanitizer.amsol import mol2amsol

# Refactored by Thua-Phong Lam, Jens Carlsson lab, Uppsala University (July, 2024)
# Based on ligand/amsol/calc_solvation.py3.csh (structure and workflow)
# Based on ligand/amsol/make_amsol71_input.py3.py (preparing input)
# Based on ligand/amsol/process_amsol_mol2.py3.py (processing AMSOL output)

"""python prepare_input.py ZINC*****.mol2 1
"""


def convert_to_ZmatMOPAC(input_file, output_file, VERBOSE=False):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mopin")
    mol = ob.OBMol()
    obConversion.ReadFile(mol, input_file)
    obConversion.WriteFile(mol, output_file)


def read_ZmatMOPAC(Zmat_file, VERBOSE=False):
    if VERBOSE:
        print("\njust entered read_ZmatMOPAC()\n")

    ZmatMOPAC_lines = {}

    with open(Zmat_file, 'r') as infile_ZmatMOPAC:
        lines = infile_ZmatMOPAC.readlines()

    # Process lines, skipping the first three
    for line_key_infile, line in enumerate(lines[3:], start=4):
        line_key_out = line_key_infile - 3

        # Modify the first three Z-matrix lines to avoid non-fatal errors in amsol7.1
        if line_key_out == 1:
            spl = line.split()
            spl[2], spl[4], spl[6] = "0", "0", "0"
        elif line_key_out == 2:
            spl = line.split()
            spl[4], spl[6] = "0", "0"
        elif line_key_out == 3:
            spl = line.split()
            spl[6] = "0"
        else:
            spl = line.split()

        # Reconstruct the line with the necessary modifications
        line = "%-2s %10.6f %2d %11.6f %2d %11.6f %2d %5d %3d %3d\n" % (
            spl[0], float(spl[1]), int(spl[2]), float(spl[3]),
            int(spl[4]), float(spl[5]), int(spl[6]),
            int(spl[7]), int(spl[8]), int(spl[9])
        )

        ZmatMOPAC_lines[line_key_out] = line

    if VERBOSE:
        print("read_ZmatMOPAC() has finished.")

    return ZmatMOPAC_lines

def create_amsol71_inputfile(output_prefix, MoleculeName, ZmatMOPAC_Data, netcharge, VERBOSE=False):

    if VERBOSE:
        print("just entered the function create_amsol71_inputfile(): ")


    # slice off the ending .ZmatMOPAC from string_Path_And_NameZmatMOPACFile by using [0:-10]:

    # open a file for the SM5.42R calculation in water solvent:
    Actual_Amsol71_InputFile_Water = open("%s" % ( output_prefix + ".in-wat") ,'w')
     
    # open a file for the SM5.42R calculation in hexadecane solvent:
    Actual_Amsol71_InputFile_Hexadecane = open("%s" % ( output_prefix + ".in-hex") ,'w')

    # the AMSOL7.1 input needs the net charge of the molecule (in its specific protonated state):
    # In earlier version, netcharge was calculated from OEChem of OpenEye. This script uses input from rdkit instead.
    if VERBOSE:
        print("netcharge of molecule in temp.mol2 (sum of partial charges):", netcharge)
    
    # write the AMSOL7.1 keywords for a SM5.42R point calculation in water to the AMSOL7.1 water input-file:
    Water_Amsol71_SM542R_Keywords = """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=WATER\n""" % netcharge
    Actual_Amsol71_InputFile_Water.write(Water_Amsol71_SM542R_Keywords)

    # write the AMSOL7.1 keywords for a SM5.42R point calculation in hexadecane to the AMSOL7.1 hexadecane input-file:
    Hexadecane_Amsol71_SM542R_Keywords = """CHARGE=%s AM1 1SCF TLIMIT=15 GEO-OK SM5.42R\n& SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93\n& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 DEV\n""" % netcharge
    Actual_Amsol71_InputFile_Hexadecane.write(Hexadecane_Amsol71_SM542R_Keywords)

    # write the name of the currently treated protonated state of the molecule into the AMSOL7.1 file
    # plus the number of atoms in the molecule
    if VERBOSE:
        print("len(ZmatMOPAC_Data) = number of atoms in molecule : ", len(ZmatMOPAC_Data))
    NumberOfAtomsInMolecule = len(ZmatMOPAC_Data)
    Molecule_Name_NrAtoms = ( "%s %d\n" % (MoleculeName, NumberOfAtomsInMolecule) )
    Actual_Amsol71_InputFile_Water.write(Molecule_Name_NrAtoms)
    Actual_Amsol71_InputFile_Hexadecane.write(Molecule_Name_NrAtoms)

    # write a blank line after the keywords block and the line showing the name of the protonated state of the molecule to the AMSOL7.1 input-files
    blank_line = "\n"
    Actual_Amsol71_InputFile_Water.write(blank_line)
    Actual_Amsol71_InputFile_Hexadecane.write(blank_line)

    # write the lines of the MOPAC Z-matrix to the AMSOL7.1 input-files
    for line_keys in ZmatMOPAC_Data:
        Actual_Amsol71_InputFile_Water.write(ZmatMOPAC_Data[line_keys])
        Actual_Amsol71_InputFile_Hexadecane.write(ZmatMOPAC_Data[line_keys])

    Actual_Amsol71_InputFile_Water.close()
    Actual_Amsol71_InputFile_Hexadecane.close()

    if VERBOSE:
        print("just finished the function create_amsol71_inputfile(). ")
 
    return


def prepare(mol2file, name, netcharge, VERBOSE=False):
    


    if VERBOSE:
        print(f"Preparing AMSOL7.1 input for {mol2file} (first: transformation to ZmatMOPAC by openbabel)")
        print("If there is any trouble, make sure that your DOCKBASE and OBABELBASE")
        print("is correctly set in ~/.cshrc or ~/.bashrc")
  

    subprocess.run(["cp", mol2file, "temp.mol2"])

    convert_to_ZmatMOPAC("temp.mol2", "temp.ZmatMOPAC", VERBOSE)
    ZmatMOPAC_data = read_ZmatMOPAC("temp.ZmatMOPAC", VERBOSE)
    create_amsol71_inputfile('temp', name, ZmatMOPAC_data, netcharge, VERBOSE)

def run(input_file, output_file, env, timeout_seconds=60, VERBOSE=False):
    AMSOLEXE = Path(__file__).parent / "amsol7.1"
    if VERBOSE: 
        print(f"Running: {AMSOLEXE} < {input_file} > {output_file}")
    try:
        subprocess.run([f"{AMSOLEXE} < {input_file} > {output_file}"], env=env, shell=True, text=True, timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        print(f"AMSOL execution timed out after {timeout_seconds} seconds")
    except subprocess.CalledProcessError as e:
        print(f"AMSOL execution failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")

def is_int(a):
    """Returns true if a can be an integer"""
    try:
        int(a)
        return True
    except:
        return False


def process_amsol_file(file, outputprefix, solvent, VERBOSE=False):
    # reads in data amsol output.
    if VERBOSE:
        print("")
        print("**** starting the function process_amsol_file() ****")
        print("")

    if not (os.path.exists(file)):
        print(file + "does not exist. \n\n Exiting script . . .")
        exit()

    # Determine if the file is gzipped and open appropriately
    splitfile = file.split('.')
    if VERBOSE:
        print(splitfile[-1])

    if splitfile[-1] == 'gz':
        input_file = gzip.open(file, 'rt')
    else:
        input_file = open(file, 'r')

    outputfilename = f"{outputprefix}{solvent}.log"
    lines = input_file.readlines()
    
    with open(outputfilename, 'w') as output:
        name = ''
        numatoms = 0
        alist = []
        total_line = []

        for line in lines:
            linesplit = line.split()  # Split on whitespace

            # Extract molecule name and atom count (eg. ZINC000007659086.1 48)
            if len(linesplit) == 2 and is_int(linesplit[1]) and not name:
                output.write(line)
                name = linesplit[0]
                numatoms = int(linesplit[1])
            # Extract the large table near the end of the AMSOL7.1 output (9 columns)
            # Extract per-atom breakdown of solvation calculation
            if len(linesplit) == 9 and is_int(linesplit[0]) and linesplit[4] != "*":
                output.write(line)
                #alist[i] contains information about the atom i
                alist.append(linesplit) 

            # Extract the first "Total: ..." line (after LS Contribution)
            elif len(linesplit) == 6 and linesplit[0] == "Total:":
                total_line = linesplit
                output.write(line)

    if VERBOSE:
        print("\n**** The function process_amsol_file() was finished. ****\n")

    return alist, total_line, name, numatoms


def diff_amsol71_files(atom_listwat,totwat,atom_listhex,tothex,name,numatom,outputprefix, VERBOSE=False):
    ## this function will compare hex with wat and write an *.solv output file:
    ## The difference "water minus hexadecane" will be calculated.
    ##
    ## The content of the *.solv file will have the following format:
    ##
    ## ZINC000000  52  0.0   -16.55   436.86     7.27    -9.28
    ## -0.1645   -0.19  19.39    0.56    0.37

    ## The first line from the above-mentioned table contains the following information :
    ## ZINC000000  52  0.0   -16.55   436.86     7.27    -9.28 means :
    ##
    ## column 1 molecule name
    ## column 2 number of atoms
    ## column 3 formal charge, i.e., net charge of the entire molecule
    ## column 4 difference of total polarization free energies for water and hex
    ## column 5 Area in (Ang**2) just from hexadecane AMSOL7.1 output-file
    ## column 6 difference of SS G-CDS values from water and hexadecane AMSOL7.1 outputs
    ## column 7 sum of column 4 and column 6: sum of the differences of the total polar and apolar atomic contributions to the solvation free energies for water and hexadecane 
    ##
    ## the following lines are the per-atom break-down:
    ## cf.: 
    ## -0.1645   -0.19  19.39    0.56    0.37 et cetera linea:
    ##
    ## column 1 is atomic partial charge just from the hexadecane (!!!) output file (hexadecane should simulate a protein-like environment: The ligand will be docked
    ##          into a protein pocket. The later use of partial charges obtained for a hexadecane environment in the electrostatic docking term is reasonable.)
    ## column 2 is the difference of the atomic polarization free energies obtained in water and hexadecane(kcal/mol): (wat - hex) 
    ## column 3 is Area (Ang**2) just from hex (is the same for water and hex. hence the source-file (whether hex or water output-files) does not really matter.)
    ## column 4 is SS G-CDS (kcal/mol) (((in amsol-mod 4 called "G-CD")))
    ## column 5 is Subtotal: (kcal/mol) (((in amsol-mod4 output called "Total Solv. free energy")))

    if VERBOSE:
        print("")
        print("**** starting the function diff_amsol71_files() ****")
        print("")

    N = len(atom_listwat)
    if len(atom_listwat) != len(atom_listhex):
        print("Error: len(atom_listwat) != len(atom_listhex):" + str(len(atom_listwat))+" != "+str(len(atom_listhex)))
        return -1

    if N != numatom:
        print("\n".join(' '.join(l) for l in atom_listwat))
        print("Error: len(atom_listwat) != numatom: " + str(N) +" != "+str(numatom))
        return -1

    fileline = ''## generate string output to write to a file
    if VERBOSE:
        print("atom_listwat[i][j] + -- + atom_listhex[i][j] +;   ::::")
        for i in range(N):
            # amsol-mod4 for j in range(0,8):
            for j in range(0,9):
               print(atom_listwat[i][j] +" -- "+ atom_listhex[i][j] +';  ', end=' ')
            print("\n", end=' ')

    # Initialize lists
    Chghex, Polhex, SAA, Apolhex, Polwat, Apolwat = ([] for _ in range(6))
    Sigmahex, diff_Pol, diff_Apol, diff_AtomicSolv = ([] for _ in range(4))

    for i in range(N):
        Chghex.append(0.0)   # list of partial atomic charges from AMSOL7.1 SM5.42R solvation calculation in hexadecane solvent
        Polhex.append(0.0)   # list of atomic polar contributions to solvation free energy: G_P obtained for hexadecane solvent
        SAA.append(0.0)      # list of atomic surface area contributions taken from hexadecane output
        Sigmahex.append(0.0) # list of sigma coefficients 
        Apolhex.append(0.0)  # list of atomic apolar contributions to solvation free energy: SS G_CDS obtained for hexadecane solvent

        Polwat.append(0.0)   # list of atomic polar contributions to solvation free energy: G_P obtained for water solvent
        Apolwat.append(0.0)  # list of atomic apolar contributions to solvation free energy: SS G_CDS obtained for water solvent

        diff_Pol.append(0.0)     # list of differences of the atomic polar contributions to solvation free energy in water and hexadecane: wat - hex
        diff_Apol.append(0.0)    # list of differences of the atomic polar contributions to solvation free energy in water and hexadecane: wat - hex 
        diff_AtomicSolv.append(0.0) # list of differences of the atomic solvation free energies (polar + apolar contributions) in water and hexadecane: wat - hex

    # Initializing the sums (over the atomic contributions):

    Chghexsum = 0.0  # sum of the partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane (!) solvent
    Polhexsum = 0.0  # sum of the atomic polar contributions to solvation free enthalpy: Polarization Free Energy (G_P) 
    SAAsum = 0.0     # sum of the atomic surface area contibutions in Angstrom
    sum_Apolhex = 0.0    # sum of the atomic apolar contributions to solvation free enthalpy: SS G_CDS (in amsol-mod4: just called G_CD, but contained S)
                   # CDS means: cavity-dispersion-solvent structure (reordering)
                   # sum_Apolhex is the same for water or hexadecane solvent. Therefore, just the hexadecane case has been considered here.
 
    FreeEnergyHex_sum = 0.0       # sum of the atomic contributions to the total solvation free energy obtained in hexadecane: In AMSOL7.1 output-table called: Subtotal

    tot_diff_Pol = 0.0            # total of the atomic polar contributions to solvation free enthalpy: wat - hex
    tot_diff_Apol = 0.0           # total of the atomic apolar contributions to solvation free enthalpy: wat - hex
    tot_diff_PolPlusApol = 0.0    # sum of polar and apolar: wat - hex



    for i in range(N):
           # alist[i][0] Atom number    (1 up to number of atoms in molecule)
           # alist[i][1] Chem. symbol
           # alist[i][2] CM2 chg.       (CM2 partial atomic charge)
           # alist[i][3] G_P (kcal/mol) (atomic polar contribution to solvation free energy)
           # alist[i][4] Area (Ang**2)  (surface area)
           # alist[i][5] Sigma (kcal/Ang**2)  (from AMSOL7.1 source code: Sigma = 1000D0*SRFACT(L)/ATAR(L))
           #                                  (i.e., Sigma[i] = 1000 * SS G_CDS[i]/Area[i], for atom i)
           #                                  ("Modeling free Energies of Solvation and Transfer", Computational Thermochemistry, Chapter 15,
           #                                  D.J. Giesen, D.G. Truhlar, 1998)
           # alist[i][6] SS G_CDS (kcal/mol) (atomic apolar contribution to solvation free energy)
           # alist[i][7] Subtotal (kcal/mol) (atomic solvation free energy (polar+apolar))
           # alist[i][8] M value


           # hexadecane:
           Chghex[i]  = float(atom_listhex[i][2])  # partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane solvent
           Polhex[i]  = float(atom_listhex[i][3])  # atomic polar contributions to solvation free enthalpy in hexadecane: Polarization Free Energy (G_P) in kcal/mol
           SAA[i]     = float(atom_listhex[i][4])  # atomic surface area contibutions in Angstrom^2
           Sigmahex[i] = float(atom_listhex[i][5]) # Sigma (kcal/Ang**2)
           # in amsol-mod4 output::: 
           #Apolhex[i] = float(atom_listhex[i][5]); # atomic apolar contributions to solvation free enthalpy in hexadecane: SS G_CDS in kcal/mol 
           Apolhex[i] = float(atom_listhex[i][6])  # atomic apolar contributions to solvation free enthalpy in hexadecane: SS G_CDS in kcal/mol 

           # water:
           Polwat[i]  = float(atom_listwat[i][3])  # atomic polar contributions to solvation free enthalpy in water: Polarization Free Energy (G_P) in kcal/mol
           # in amsol-mod4 output::: 
           #Apolwat[i] = float(atom_listwat[i][5]) # atomic apolar contributions to solvation free enthalpy in water: SS G_CDS in kcal/mol
           Apolwat[i] = float(atom_listwat[i][6])  # atomic apolar contributions to solvation free enthalpy in water: SS G_CDS in kcal/mol

           ## sums over atomic contributions:

           Chghexsum += Chghex[i] # sum of the partial charges (CM2: Truhlar's charge model 2 charges) of each atom in hexadecane (!) solvent
           Polhexsum += Polhex[i] # sum of the atomic polar contributions to solvation free enthalpy: Polarization Free Energy (G_P)
           SAAsum += SAA[i]          # sum of the atomic surface area contibutions in Angstrom
                                     # is independent of the solvent under consideration
           # !!! former apolsum
           sum_Apolhex += Apolhex[i]    # sum of the atomic apolar contributions to solvation free enthalpy: SS G_CDS from hexadecane calculation (!!!!)
           # !!! former energy_sum
           # in amsol-mod4 output::: 
           #FreeEnergyHex_sum = FreeEnergyHex_sum + float(atom_listhex[i][6]); # sum of the atomic contributions to the total solvation free energy
                                                                               # for hexadecane (!!!) solvent: In AMSOL7.1 output-table called: Subtotal (in kcal/mol)
           FreeEnergyHex_sum += float(atom_listhex[i][7]); # sum of the atomic contributions to 
                                                           # for hexadecane (!!!) solvent: In AMSOL7.1 output-table called: Subtotal (in kcal/mol)
    
    # The information written in the "Total: ..." line close to the end of the hexadecane AMSOL7.1 output:
    #
    # tothex[0] = "Total:"
    # tothex[1] = CM2 chg.
    # tothex[2] = G_P, i.e. polar contribution to solvation free energy or also called Polarization free energy (in kcal/mol)
    # tothex[3] = Area (Ang**2) (in amsol-mod4: Area (CD) (Ang**2)): sum of all atomic contributions (!)
    # tothex[4] = SS G_CDS in kcal/mol (in amsol-mod4: G-CD)
    # tothex[5] = Subtotal in kcal/mol, i.e. sum of all atomic polar (G_P[i]) and apolar (SS G_CDS[i]) contributions to the solvation free energy
    #            Subtotal = solvation free energy
    #            Subtotal = "(5)  G-P-CDS(sol) = G-P(sol) + G-CDS(sol) = (2) + (4)             -XX.XXX kcal" 

    cs_coeff = (float(tothex[4]) - sum_Apolhex)/float(tothex[3])   # cs_coeff = total LS contribution / total AreaSAA
    #           (total SS G_CDS - sum of atomic all atomic apolar contributions to free energy of solvation)/ total Area
    #           = total CS/LS contribution divided by total Area


    sum_csTimesSAA = 0.0

    # !!! WATER minus HEXADECANE: performing wat-hex difference !!! : 
    for i in range(N): # do substraction
           diff_Pol[i] = Polwat[i] - Polhex[i]                               # diff_Pol[i] is the difference of the atomic polar contribution to solvation free energy                                                                            # obtained in water vs. the one obtained in hexadecane
                                                                             
           tot_diff_Pol = tot_diff_Pol + diff_Pol[i]                         # tot_diff_Pol is the sum of all differences of atomic polar contributions to solvation
                                                                             # free energy obtained in water vs. in hexadecane

           sum_csTimesSAA = sum_csTimesSAA + cs_coeff * SAA[i]               # sum_i cs_coeff * SAA[i] = sum_i (LS contribution/ total Area) * atomic Surface Area SAA[i] = LS contribution

                                                                             # - (Apolhex[i] + cs_coeff * SAA[i]): minus since: wat MINUS hex. cs_coeff and LS contribution only present in hex (Large Solvent)
           diff_Apol[i] = Apolwat[i] - (Apolhex[i] + cs_coeff * SAA[i])      # diff_Apol[i] is the difference of the atomic apolar contribution to solvation free energy
                                                                             # obtained in water vs. the one obtained in hexadecane minus cs_coeff * SAA[i]  
                                                                             
           tot_diff_Apol = tot_diff_Apol + diff_Apol[i]                      # tot_diff_Apol is sum of all atomic diff_Apol[i]
                                                                             
           diff_AtomicSolv[i] = diff_Pol[i] + diff_Apol[i]                   # diff_AtomicSolv[i] is the sum of the atomic diff_Pol[i] and diff_Apol[i] (see above for diff_Pol[i] and diff_Apol[i])
                                                                             # diff_AtomicSolv[i] considers the difference between water and hexadecane results for each atom!!!

           tot_diff_PolPlusApol = tot_diff_PolPlusApol + diff_AtomicSolv[i]  # tot_diff_PolPlusApol is the sum over all atomic contributions to the difference of solvation free energies
                                                                             # in water and hexadecane, i.e. the difference
                                                                             # of the final solvation free energies obtained in water and hexadecane 

    ## write out the solvation file 
    file = open(outputprefix+'.solv','w')
    file.write( "%s %3d %4.1f %8.2f %8.2f %8.2f %8.2f\n" % (name,numatom,float(totwat[1]),tot_diff_Pol,float(totwat[3]),tot_diff_Apol,tot_diff_PolPlusApol))  # formal charge is the same for water or hexadecane
    for i in range(N):
           file.write("%8.4f%8.2f%7.2f%8.2f%8.2f\n" % (Chghex[i],diff_Pol[i],SAA[i],diff_Apol[i],diff_AtomicSolv[i]))
    file.close()

    if VERBOSE:
        print("\n**** The function diff_amsol71_files() was finished. ****\n")

    return 0

#################################################################################################################
#################################################################################################################

def modify_charges_mol2_file(mol2file, atom_list_hex, outputprefix, VERBOSE=False):
    ## read in mol2 file

    if VERBOSE:
        print("")
        print("**** starting the function modify_charges_mol2_file() *****")
        print("")
        print("     CM2 charges from an AMSOL7.1 SM5.42R calculation in hexadecane (!!!) ")
        print("     are written to mol2-file. Former charges are overwritten.")

    mol = mol2amsol.read_Mol2_file(mol2file)[0] 

    n = len(atom_list_hex)
    if n != len(mol.atom_list):
       print("Error: n != len(mol.atom_list) : " + str(n) + " !=" + str(len(mol.atom_list)))
       return -1

    
    for i in range(n):
        charge = atom_list_hex[i][2]
        if VERBOSE:
            print(mol.atom_list[i].Q, charge)
        mol.atom_list[i].Q = float(charge)

    filename = outputprefix + '.mol2'
    mol2amsol.write_mol2(mol,filename)
    
    if VERBOSE:
        print("")
        print("**** The function modify_charges_mol2_file() was finished. ****")
        print("")

    return 0
    
def process_output(wat_file, hex_file, mol2file, output_prefix, VERBOSE=False):
    atom_list_wat,tot_wat,name_wat,numat_wat = process_amsol_file(wat_file,output_prefix,"wat")
    atom_list_hex,tot_hex,name_hex,numat_hex = process_amsol_file(hex_file,output_prefix,"hex", VERBOSE = VERBOSE)
    

    # tot_wat and tot_hex are lists:
    # ( see def process_amsol_file(...) above)
    #
    # tot_wat[0] or tot_hex[0] : total_line[0] = "Total:"
    # tot_wat[1] or tot_hex[1] : total_line[1] = CM2 chg.
    # tot_wat[2] or tot_hex[2] : total_line[2] = G_P, i.e. polar contribution to solvation free energy or also called Polarization free energy (in kcal/mol)
    # tot_wat[3] or tot_hex[3] : total_line[3] = Area (Ang**2) (in amsol-mod4: Area (CD) (Ang**2))
    # tot_wat[4] or tot_hex[4] : total_line[4] = SS G_CDS in kcal/mol (in amsol-mod4: G-CD)
    # tot_wat[5] or tot_hex[5] : total_line[5] = Subtotal in kcal/mol, i.e. sum of all atomic polar (G_P[i]) and apolar (SS G_CDS[i]) contributions to the solvation free energy
    #                            Subtotal = solvation free energy
    #                            Subtotal = "(5)  G-P-CDS(sol) = G-P(sol) + G-CDS(sol) = (2) + (4)             -XX.XXX kcal"
    

    error_signal = 0 # 0 means no error -1 means error occured
    if (name_hex != name_wat or numat_hex != numat_wat):
        print("Error: Name or Atom counts do not agree")
    elif VERBOSE:
        print("Wat-name      = " + wat_file) 
        print("Wat-atom cout = " + str(numat_wat))
        print("Hex-name      = " + hex_file) 
        print("Hex-atom cout = " + str(numat_hex))

    if VERBOSE:
        print("")
        print("just before diff_amsol71_files() function")
        print("")
    error_signal = diff_amsol71_files(atom_list_wat, tot_wat, atom_list_hex, tot_hex, name_wat, numat_wat, output_prefix)

    if VERBOSE:
        print("")
        print("just before modify_charges_mol2_file() function")
        print("")
    if error_signal == 0: modify_charges_mol2_file(mol2file, atom_list_hex, output_prefix) 

    if VERBOSE:
        print("")
        print("**** The main program in process_amsol71_mol2.py was finished for ****")
        #print("%s and %s." % (filenamewat, filenamehex))
        print("*******************************************************************")
        print("")
    #
    return error_signal