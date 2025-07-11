#!/usr/bin/env python3.7

#Ryan G. Coleman
#reads in file containing clash definitions. writes defaults.

#import string
from __future__ import print_function
import sys
from msani.db2.geometry import distL2Squared3
from collections import defaultdict
import operator

class Clash(object):
  '''holds parameters that determine what a clashed conformation contains.
  rules have the format ('min|max', bond, cmp, "X", "Y", dist) which is:
  min|max - whether the constraint is a minimum or maximum
  bond - how many bonds are between the atoms.
  cmp - how to use the bond distance. -1 means there must be less than that many
    0 means there must be exactly that many and
    1 means there must be more than that many bonds between the atoms
  X - atom type X, * means all
  Y - atom type Y, * means all
  dist - float that is the distance constraint'''
  #old rules, necessary for mix'n'match, but not for simple hydroxyl rotations
  #  rulesDefault = [("max", 1, 0, "*", "*", 2.2),
  #                  ("min", 3, 1, "O", "O", 2.0),
  #                  ("min", 3, 1, "*", "*", 1.50),
  #                  ("min", 1, 0, "O", "O", 2.0),
  #                  ("min", 1, 0, "*", "*", 0.95)]
  #  rulesDefault = [("min", 2, 1, "*", "*", 1.70)]
  rulesDefault = [("min", 2, 1, "H", "H", 1.60)]   # only H-H!!

  def __init__(self, parameterFileName=None):
    '''constructs from defaults or reads from file'''
    if parameterFileName is not None:
      parameterFile = open(parameterFileName, 'r')
      self.rules = []
      try:
        for line in parameterFile:
          #tokens = string.split(line)
          tokens = line.split()
          self.rules.append((
              tokens[0], int(tokens[1]), int(tokens[2]),
              tokens[3], tokens[4], float(tokens[5]), float(tokens[5])**2.))
      except StopIteration:
        pass  # EOF
    else:  # no parameter file, use defaults
      self.rules = []
      for defaultRule in self.rulesDefault:
        ruleWithSquared = list(defaultRule)
        ruleWithSquared.append(defaultRule[5]**2.)
        self.rules.append(tuple(ruleWithSquared))

  def printParameters(self):
    '''prints to standard out the parameters used in a readable format'''
    for rule in self.rules:
      for part in rule[0:6]:  # don't print out the distance squared term
        print(part, end=' ')
      print("")  # force newline

  def decideDistanceRules(self, mol2data, xyzData):
    natoms = len(xyzData)

    for rule in self.rules:
      const = rule[0]
      bondc = rule[1]
      compr = operator.lt if rule[2] == -1 else operator.eq if not rule[2] else operator.gt
      distc = rule[5]
      distc2 = distc**2
      typeA = rule[3]
      typeB = rule[4]

      if typeA == typeB:

        # ye olde fashioned way
        atoms = [atom for atom in range(natoms) if 0 == mol2data.atomType[atom].find(typeA) or typeA == "*"]
        ntypeatoms = len(atoms)

        for atomi in range(ntypeatoms):
          for atomj in range(atomi+1, ntypeatoms):
            atomA = atoms[atomi]
            atomB = atoms[atomj]
            if not compr(mol2data.bondsBetweenActual(atomA, atomB), bondc):
              continue
            dist = distL2Squared3(xyzData[atomA], xyzData[atomB])
            if ((const == "min" and dist < distc2) or (const == "max" and dist > distc2)):
              return True

        """
        # was trying to figure out a way to do this with the bucketing algorithm
        # unfortunately I couldn't think of a great way for this to work out, the clash deciding problem is a different beast than point clustering
        # I've settled for optimizing the current algorithm the best I can. Since we only work with one rule, H-H, it is faster to skip calculating distances between all atom pairs
        # and just calculate distances between the hydrogen atoms. This seems to speed up clash deciding significantly. If we ever go back to using more than one rule, this algorithm
        # will have to be reworked again
        bucketer = buckets2.buckets2(distc)
        confClusters = {}
        setToConfs = {}

        atomsxyz = [xyzData[atom] for atom in range(natoms) if 0 == mol2data.atomType[atom].find(typeA) or typeA == "*"]
        ntypeatoms = len(atomsxyz)
        
        npos, confNum = bucketer.bucket(atomsxyz, confClusters=confClusters, setToConfs=setToConfs)

        if const == "min" and npos != ntypeatoms:
          return True
        
        if const == "max" and npos != 1:
          return True
        """
          
      # otherwise we need to do things the old fashioned way, checking distance between pairs
      else:

        if (typeA == "*" and typeB != "*") or (typeA != "*" and typeB == "*"):
          typeDef = typeB if typeA == "*" else typeA

          atomsDef = [atom for atom in range(natoms) if 0 == mol2data.atomType[atom].find(typeDef)]
          ndefatoms = len(atomsDef)
          natomsminusdef = natoms - ndefatoms
          atomsAllMinusDef = [atom for atom in range(natoms) if 0 != mol2data.atomType[atom].find(typeDef)]

          for atomi in range(ndefatoms):
            for atomj in range(atomi+1, ndefatoms):
              atomA = atomsDef[atomi]
              atomB = atomsDef[atomj]
              if not compr(mol2data.bondsBetweenActual(atomA, atomB), bondc):
                continue
              dist = distL2Squared3(xyzData[atomA], xyzData[atomB])
              if (const == "min" and dist < distc2) or (const == "max" and dist > distc2):
                return True

          for atomA in atomsDef:
            for atomB in atomsAllMinusDef:
              if not compr(mol2data.bondsBetweenActual(atomA, atomB), bondc):
                continue
              dist = distL2Squared3(xyzData[atomA], xyzData[atomB])
              if (const == "min" and dist < distc2) or (const == "max" and dist > distc2):
                return True

        else:

          atomsA = [atom for atom in range(natoms) if 0 == mol2data.atomType[atom].find(typeA) or typeA == "*"]
          atomsB = [atom for atom in range(natoms) if 0 == mol2data.atomType[atom].find(typeA) or typeA == "*"]

          for atomA in atomsA:
            for atomB in atomsB:
              if not compr(mol2data.bondsBetweenActual(atomA, atomB), bondc):
                continue
              dist = distL2Squared3(xyzData[atomA], xyzData[atomB])
              if ((const == "min" and dist < distc2) or (const == "max" and dist > distc2)):
                return True

    return False

  # old clash decider algorithm
  # here for reference's sake
  def decide(self, mol2data, xyzData):
    '''mol2data is the mol2.Mol2 object. xyzData is a list of coords.
    use self.rules to return True (clashed) or False (not clashed).'''
    dists = defaultdict(list)  # format is atomNum -> (otherNum, dist, bondDist)
               #all dists in list are euclidean distance squared
    atomNums = list(range(len(xyzData)))
    atomNums.sort()
    for atomNumOne in atomNums:
      for atomNumTwo in atomNums:
        if atomNumTwo > atomNumOne:
          thisDist = distL2Squared3(xyzData[atomNumOne], xyzData[atomNumTwo])
          bondDist = mol2data.bondsBetweenActual(atomNumOne, atomNumTwo)
          dists[atomNumOne].append((atomNumTwo, thisDist, bondDist))
          dists[atomNumTwo].append((atomNumOne, thisDist, bondDist))
    for rule in self.rules:
      #match atom types first
      for atomNum in atomNums:
        if rule[3] == "*" or \
            0 == mol2data.atomType[atomNum].find(rule[3]):
            # fix for py3-3.7
            #0 == string.find(mol2data.atomType[atomNum], rule[3]):
          for dist in dists[atomNum]:  # for every distance
            if rule[4] == "*" or \
                0 == mol2data.atomType[dist[0]].find(rule[4]):
                # fix for py3-3.7
                #0 == string.find(mol2data.atomType[dist[0]], rule[4]):
              brokeRule = False
              if rule[0] == "max":  # is a max distance constraint
                if dist[1] > rule[6]:  # broke the rule
                  brokeRule = True
              elif rule[0] == "min":  # is a min distance constraint
                if dist[1] < rule[6]:  # broke the rule
                  brokeRule = True
              if brokeRule:  # check to make sure actually broken
                cmp = (dist[2] > rule[1]) - (dist[2] < rule[1])
                #if not cmp(dist[2], rule[1]) == rule[2]:  # this amounts
                if not cmp == rule[2]:  # this amounts
                    #to checking to see if the right number of bonds lie between
                    #the atoms in question.
                  brokeRule = False
                if brokeRule:  # rule has been broken so there is a clash
                  #print rule, atomNum, dist #debug the rules broken
                  return True  # can quit after first broken rule
    #if everything passed, return False indicating no clashes
    return False

# fix for py3-3.7
#if -1 != string.find(sys.argv[0], "clash.py"):
if -1 != sys.argv[0].find("clash.py"):
  #if program is called from the command line, assume user wants a copy of the
  #default parameter file written to standard out. this is the only command use.
  #usually this will be imported and run from somewhere else.
  if len(sys.argv) > 1:
    Clash(sys.argv[1]).printParameters()
  else:
    Clash().printParameters()
