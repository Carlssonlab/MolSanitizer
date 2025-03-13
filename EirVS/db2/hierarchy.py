#!/usr/bin/env python3.7

#Ryan G. Coleman
#uses mol2 file to generate a hierarchy

#import string
import sys
from .unionfind2 import unionFind
from .geometry import distL2Squared3
from . import geometry
from . import buckets
import gzip
import operator
import math
import time
from . import shortestpaths
import itertools
from . import buckets2
import io

def printClusterHelper(clusterList):
  '''stupid function used for debugging, prints list of pymol out.???.mol2 lines
  to copy/paste and run to see what the clusters are.
  run mol2hydroxyls.py -r and mol2tomultimol2.py first to get out.???.mol2 files
  '''
  for clusters in clusterList:
    print("pymol ", end=" ")
    for conf in clusters:
      #print("out." + string.zfill(conf, 3) + ".mol2 ", end=' ')
      print("out." + conf.zfill(3) + ".mol2 ", end=" ")
    print(" ")

def computeBreaks(limitError, options):
  '''3 diff requirements, make sure we break it into enough pieces to meet them
  all.'''
  #have to break the atomXyz into multiple sets so the hierarchy isn't too big
  try:
    breaksS = int(math.ceil(limitError.getSets() / float(options.limitset)))
  except TypeError:  # means None was used
    breaksS = 1
  try:
    breaksC = int(math.ceil(limitError.getConfs() / float(options.limitconf)))
  except TypeError:  # means None was used
    breaksC = 1
  try:
    breaksX = int(math.ceil(limitError.getCoords() / float(options.limitcoord)))
  except TypeError:  # means None was used
    breaksX = 1
  #print breaksS, breaksC, breaksX  # see which breaks is higher
  breaks = max(breaksS, breaksC, breaksX)  # use the max of any of these
  return breaks

class TooBigError(Exception):
  '''error raised when the hierarchy has too many conformations of input
  after the hydroxyls have been rotated.'''

  def __init__(self, confs, sets, coords):
    self.confs = confs
    self.coords = coords
    self.sets = sets

  def __str__(self):
    return repr(self.confs) + ", " + repr(self.coords) + ", " + repr(self.sets)

  def getConfs(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.confs

  def getCoords(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.coords

  def getSets(self):
    '''actually used to figure out how many sub-groups to split input confs'''
    return self.sets

class Hierarchy(object):
  '''uses data from a mol2 file to make a hierarchy of conformations.
  the following constants are used when writing out the confs/groups and are
  based on the 80 character limit in fortran. yeah seriously.
  they might change if something serious happens but it is better that they
  are here than hardcoded several times later
  these are floats so that the division works'''
  grGrPerLine = 17.  # group -> group children per line in output
  grCoPerLine = 9.  # group -> conf
  coCoPerLine = 9.  # conf -> conf
  coSePerLine = 8.  # conf -> set

  def __init__(
      self, mol2data, clashDecider, tolerance=0.001, verbose=False,
      timeit=False, limitset=9999999999, limitconf=9999999999,
      limitcoord=9999999999, solvdata=None):
    '''takes a mol2data class as input. makes a hierarchy.'''
    if solvdata is not None:
      self.solvdata = solvdata
    if timeit:
      startTime = time.time()
    #first step is to count the number of positions each atom has.
    #the tolerance is taken into account here and only here.
    totalCoords = len(mol2data.atomXyz) * len(mol2data.atomXyz[0])
    if verbose:
      print("total number of sets (complete confs):", len(mol2data.atomXyz))
    if len(mol2data.atomXyz) > limitset:  # quit now, way too many sets
      if verbose:
          print("Total number of sets too high %d > %d" % (len(mol2data.atomXyz), limitset))
      raise TooBigError(None, len(mol2data.atomXyz), totalCoords)

    self._getRigidStructures(len(mol2data.atomXyz[0]), mol2data.atomBonds, verbose)
    if timeit:
      structureTime = time.time()
      print("time to get rigid structures:", structureTime-startTime)

    #if len(mol2data.atomXyz) > 50:
    if verbose:
      print("using faster count positions algorithm for large data")
    self._countPositions(mol2data.atomXyz, tolerance, verbose)
    #else:
    #  if verbose:
    #    print("using default count positions algorithm for smaller data")
    #  self._countPositionsFewPoints(mol2data.atomXyz, tolerance)
    if timeit:
      countTime = time.time()
      print("time to count unique positions:", countTime-structureTime)
    #if verbose:
    #  print("unique positions, atoms:", self.posCount, len(mol2data.atomXyz))
    if totalCoords > limitcoord:
      if verbose:
          print("Total number of coords too high %d > %d" % (totalCoords, limitcoord))
      raise TooBigError(None, len(mol2data.atomXyz), totalCoords)
      #this breaks out of the init stage, needs fewer confs to be passed in.
    #the rigid component is the biggest set of bonded non-moving atoms
    self._findRigidComponent(mol2data.atomBonds)  # also uses self.posCount
    #if timeit:
    #  rigidTime = time.time()
    #  print("time to find rigid component:", rigidTime-countTime)
    #if verbose:
    #  print("rigid atoms, others:", self.rigidComponent, self.atomsNotAssigned)
    #new algorithm, find bonded atoms that move together, put in conformations
    self._findRigidHeavy(mol2data.atomType)
    self.heavyAtomNums = None
    self._setHeavy(mol2data.atomType)
    # self._findConformations(mol2data.atomBonds, mol2data.atomXyz) # not needed anymore, moved conformations, set finding, and position assigning to _countPositions
    # self._findSets()  # puts conformations in sets
    #if timeit:
    #  flexTime = time.time()
    #  print("time to find flexible components:", flexTime-countTime)
    if verbose:
      print("total number of confs:", self.numConfs)
    if self.numConfs > limitconf:
      raise TooBigError(
          self.numConfs, len(mol2data.atomXyz), totalCoords)
      #this breaks out of the init stage, needs fewer confs to be passed in.
    #now want to actually put atom positions into hierarchy groups
    # self._assignCoords(mol2data.atomXyz) # moved into _countPositions
    #if timeit:
    #  assignCoordsTime = time.time()
    #  print("time to assign coords:", assignCoordsTime-flexTime)
    self._identifyClashSetnums(clashDecider, mol2data)
    if timeit:
      afterClash = time.time()
      print("time to identify clash sets:", afterClash-countTime)
    if verbose:
      print("number of broken/clashed sets:", len(self.brokenSets))
    #the mol2data is needed during output so save it.
    self.mol2data = mol2data
    #if timeit:
    #  afterXyz = time.time()
    #  print("time to identify conf atoms:", afterXyz - flexTime)
    self.clusters = None  # used to detect if clustering/clouding was done
    #self._makeClouds()  # highest level of ligand sampling
    #if timeit:
    #  afterClouds = time.time()
    #  print("time to make clouds:", afterClouds - afterClash)

    if timeit:
      endTime = time.time()
      print("time spent processing hierarchy total", endTime - startTime)

  # Benjamin Tingle 2/15/2021
  # New function, _getRigidStructures
  # gets all independently moving structures within the atom
  # creates a map from atom -> structure
  # independently moving structures contain atoms that fullfill one of the following conditions:
  #   a. connected via rigid bonds (i.e not single bonds)
  #     1. if an atom only has one bond neighbor, that bond is considered rigid, even if it is a single bond
  #   b. within a graph cycle together (graph cycles within mols are considered rigid i.e ring systems)
  # each conf should only contain atoms that move together, so we generate this information to assist in conf creation
  def _getRigidStructures(self, natoms, atomBonds, verbose):

    cycleCount = 0
    cycles = [[] for i in range(natoms)]
    parent = [0 for i in range(natoms)]
    visited = [0 for i in range(natoms)]

    def findCycles(atom, prev, parent, visited, cycles):
      nonlocal cycleCount
      if visited[atom] == 2:
        return
      if visited[atom] == 1: # found cycle!
        curr = prev
        while curr != atom:
          cycles[curr].append(cycleCount)
          curr = parent[curr]
        cycles[curr].append(cycleCount)
        cycleCount += 1
        return
      visited[atom] = 1
      parent[atom] = prev
      for otheratom, bondtype in atomBonds[atom]:
        if otheratom == prev:
          continue
        findCycles(otheratom, atom, parent, visited, cycles)
      visited[atom] = 2

    # finds all cycles in the atom e.g ring systems
    findCycles(0, -1, parent, visited, cycles)

    # adds rigid atom bonds as length 2 cycles to the cycles graph
    for atom in range(natoms):
      if len(atomBonds[atom]) == 1:
        cycles[atom].append(cycleCount)
        cycles[atomBonds[atom][0][0]].append(cycleCount)
        cycleCount += 1
        continue
      for otheratom, bondtype in atomBonds[atom]:
        if otheratom < atom:
          continue
        skip = False
        other = cycles[otheratom]
        for cycle in cycles[atom]:
          if cycle in other:
            skip = True
            break
        if not skip and bondtype != "1":
          cycles[atom].append(cycleCount)
          cycles[otheratom].append(cycleCount)
          cycleCount += 1
    
    intersections = [set() for i in range(cycleCount)]

    # find all points of intersection for each cycle
    for atomCycles in cycles:
      lcycle = len(atomCycles)
      for i in range(lcycle):
        for j in range(i+1, lcycle):
          intersections[atomCycles[i]].add(atomCycles[j])
          intersections[atomCycles[j]].add(atomCycles[i])

    rigidMap = [0 for i in range(cycleCount)]
    visited = [False for i in range(cycleCount)]
    rigidStructureCount = 0

    def findRigidStructures(cycle, prev, intersections, rigidMap, visited):
      nonlocal rigidStructureCount
      if visited[cycle]:
        return
      visited[cycle] = True
      rigidMap[cycle] = rigidStructureCount
      for intersect in intersections[cycle]:
        if intersect == prev:
          continue
        findRigidStructures(intersect, cycle, intersections, rigidMap, visited)
    
    # merge intersecting cycles into a rigid structures map (cycle -> structure)
    for cycle in range(cycleCount):
      if not visited[cycle]:
        findRigidStructures(cycle, -1, intersections, rigidMap, visited)
        rigidStructureCount += 1

    self.rigidStructures = [0 for i in range(natoms)]

    # resolve the cycle -> structure map into an atom -> structure map
    for atom in range(natoms):
      if len(cycles[atom]) > 0:
        self.rigidStructures[atom] = rigidMap[cycles[atom][0]]
      else:
        # if the atom is not part of any cycles (example: single atom connected to other structures by two single bonds)
        # then it is a rigid structure itself
        self.rigidStructures[atom] = rigidStructureCount
        rigidStructureCount += 1

    self.rigidStructureCount = rigidStructureCount
    if verbose:
      print("RIGID_STRUCTURES", self.rigidStructures)

  def _countPositions(self, xyzData, tolerance, verbose=False):
    '''
    for a list of list of xyz data, count the number of positions each
    atom takes based on the tolerance and the distance. tolerance is compared
    to the euclidean difference squared to determine if a position is equal.
    actually uses a clustering algorithm and uses a unionfind data structure.
    
    2/11/2021, Benjamin Tingle
    Reworked _countPositions massively. Uses a new bucketing algorithm that doesn't require looping through all possible buckets, also removed the unionfind thing
    Uses hashed buckets instead of an enormous array of buckets
    Also reworked it so that sets, confs, atom output coordinates are resolved in the _countPositions function rather than outside
    Confs are split along predefined rigid segments in the molecule, much easier to handle than attempting to figure out atom connectivity for each conf
    '''

    nmol2s = len(xyzData)
    natoms = len(xyzData[0])

    confClusters = {}
    posCount = []
    posTotal = 0
    confNum = 0

    bucketer = buckets2.buckets2(tolerance)

    for atom in range(natoms):

      xyzData_t = [xyzData[mol2][atom] for mol2 in range(nmol2s)]

      # check out buckets2.py for the new clustering algorithm
      npos, confNum = bucketer.bucket(xyzData_t, atomId=atom, confNum=confNum, confClusters=confClusters)
      posTotal += npos

      posCount.append(npos)

    self.confAtoms = []
    self.confInput = []
    self.outAtomOrigAtom = [0 for i in range(posTotal)]
    self.outAtomConfNum  = [0 for i in range(posTotal)]
    self.outAtomXYZ      = [0 for i in range(posTotal)]
    self.confNumAtomList = []
    self.setToConfs = {i : [] for i in range(nmol2s)}

    globalAtomCnt = 0
    confNum_act = 0
    for tupleInput, confInfo in sorted(confClusters.items(), key=lambda x:-len(x[0])): # sort confclusters such that the rigid conf is evaluated first
      confNum_t, atoms, xyzlist = confInfo

      rs_prev = None
      # split confs if they contain atoms in different rigid structures i.e the atoms separated by a rotatable bond
      # the previous method explored atom bonds within a full cluster match to find which pieces are connected and which aren't, the purpose being to merge atoms that "move together" into their own conf(s)
      # problem is, this allowed atoms that don't actually "move together" to nonetheless be assigned the same conf because rotatable bonds weren't taken into account
      # I am unsure of the purpose of having atoms that move together in the same conf, but I will not question it, just try and improve it
      rsatoms = [(a, self.rigidStructures[a]) for a in atoms]
      for i, (atom, rs) in enumerate(sorted(rsatoms, key=lambda x:x[1])):
        
        if rs_prev == rs:
          self.confAtoms[confNum_act-1].append(atom)

        elif rs_prev != rs:
          self.confAtoms.append([atom])
          self.confInput.append(tupleInput)

          if rs_prev != None:
            self.confNumAtomList[confNum_act-1][1] = globalAtomCnt - 1
          self.confNumAtomList.append([globalAtomCnt, 0])

          for setno in tupleInput:
            self.setToConfs[setno].append(confNum_act)

          confNum_act += 1

        self.outAtomOrigAtom[globalAtomCnt] = atoms[i]
        self.outAtomConfNum[globalAtomCnt] = confNum_act
        self.outAtomXYZ[globalAtomCnt] = xyzlist[i]

        globalAtomCnt += 1
        rs_prev = rs
      self.confNumAtomList[confNum_act-1][1] = globalAtomCnt - 1

    self.outAtoms = globalAtomCnt
    self.numConfs = confNum_act
    self.posCount = posCount

    if verbose:
      print("POS_COUNT", self.posCount)

  def _countPositionsFewPoints(self, xyzData, tolerance):
    '''for a list of list of xyz data, count the number of positions each
    atom takes based on the tolerance and the distance. tolerance is compared
    to the euclidean difference squared to determine if a position is equal.
    actually uses a clustering algorithm and uses a unionfind data structure.'''
    self.posCount = []
    self.posClusters = []  # just save all the data since we made it
    self.posClusterLists = []  # just save all the data since we made it
    tolerance2 = tolerance ** 2.  # square the tolerance since it is compared
    for oneSet in range(len(xyzData[0])):  # goes from 0 to atom count
      clusters = unionFind()
      xyzList = []
      for oneIndex in range(len(xyzData)):  # 0 to number of positions (mol2#s)
        clusters.find(oneIndex)  # initiate each position
        xyzList.append(xyzData[oneIndex][oneSet])
      for oneIndex in range(len(xyzData)):  # 0 to positions
        oneXyz = xyzList[oneIndex]
        for twoIndex in range(oneIndex+1, len(xyzData)):
          # count from oneIndex to positions
          if geometry.distL2Squared3(oneXyz, xyzList[twoIndex]) < tolerance2:
            clusters.union(oneIndex, twoIndex)
      tempLists = clusters.toLists()
      self.posCount.append(len(tempLists))
      self.posClusters.append(clusters)
      self.posClusterLists.append(tempLists)

  def _findRigidComponent(self, atomBonds):
    '''uses bond and position count information to find largest set of atoms
    that don't move. this is the rigid component. set into self.rigidComponent
    also find the complement of atomnums and the rigid component and set into
    self.atomsNotAssigned for use later'''
    clusters = unionFind()
    for atomNum in range(len(self.posCount)):
      if 1 == self.posCount[atomNum]:
        for otherNum, bondType in atomBonds[atomNum]:
          if 1 == self.posCount[otherNum]:
            clusters.union(atomNum, otherNum)
    maxSize = 0
    maxCluster = None
    clusterLists = clusters.toLists()
    for clusterList in clusterLists:
      if len(clusterList) > maxSize:
        maxSize = len(clusterList)
        maxCluster = clusterList
    self.rigidComponent = maxCluster
    self.atomsAssigned = set(self.rigidComponent)
    self.atomsNotAssigned = set()
    for atomNum in range(len(self.posCount)):
      if atomNum not in self.rigidComponent:
        self.atomsNotAssigned.add(atomNum)

  def _findRigidHeavy(self, atomTypes):
    '''counts the heavy atoms in the rigid component and puts in
    self.heavyRigidCount'''
    self.heavyRigidCount = 0
    self.heavyRigidAtomNums = []
    for atomNum in self.atomsAssigned:
      if atomTypes[atomNum].find('H') == -1:
        self.heavyRigidAtomNums.append(atomNum)
        self.heavyRigidCount += 1
    #print self.heavyRigidCount

  def _setHeavy(self, atomTypes):
    '''for all atoms, finds the heavy ones, put in self.heavyAtomNums, return'''
    if self.heavyAtomNums is None:  # only do this once, it never changes
      self.heavyAtomNums = []
      for atomNum in range(len(atomTypes)):
        if atomTypes[atomNum].find('H') == -1:
          self.heavyAtomNums.append(atomNum)
    return self.heavyAtomNums

  ### DEPRECATED
  def _findConformations(self, atomBonds, xyzData):
    '''uses bond and xyzs to figure out what sets of neighboring atoms move
    together and assign them to conformations and assign each set a specific
    bunch of conformations.
    self.rigidComponent is the list of atom numbers for the rigid comp
    self.atomsAssigned is the set of atom numbers for the rigid comp (@start)
    self.atomsNotAssigned is the rest of the atom numbers'''
    self.confNums = [1]  # rigid starts
    self.confAtoms = {}  # maps to atom numbers
    self.confAtoms[1] = list(self.atomsAssigned)
    self.confInput = {}  # maps to the input xyz lists
    self.confInput[1] = list(range(len(xyzData)))
    confClusters = {}
    for atomNum in self.atomsNotAssigned:
      for listInputs in self.posClusterLists[atomNum]:
        tupleInputs = tuple(listInputs)  # can't use lists as keys
        if not confClusters.get(tupleInputs):
          confClusters[tupleInputs] = unionFind()
        confClusters[tupleInputs].find(atomNum)  # in case of singletons
        for otherNum, bondType in atomBonds[atomNum]:
          if listInputs in self.posClusterLists[otherNum]:
            confClusters[tupleInputs].union(atomNum, otherNum)
    for tupleInputs, clusters in confClusters.items():
      for atomLists in clusters.toLists():
        #make a conf for each
        thisConf = self.confNums[-1] + 1
        self.confAtoms[thisConf] = atomLists
        self.confInput[thisConf] = tupleInputs
        self.confNums.append(thisConf)
    #print self.confNums, self.confAtoms, self.confInput
    #that's it, confs have been built

  ### DEPRECATED
  def _findSets(self):
    '''puts conformations together into sets'''
    self.setToConfs = {}  # maps set numbers to conf lists
    for confNum in self.confNums:
      for tupleInput in self.confInput[confNum]:
        if tupleInput not in self.setToConfs:
          self.setToConfs[tupleInput] = []
        self.setToConfs[tupleInput].append(confNum)
    #print self.setToConfs
    #self.setToConfs contains relevant mapping

  ### DEPRECATED
  def _assignCoords(self, xyzData):
    '''for each conf (including rigid) find atom positions for each atom'''
    self.outAtoms = 0  # counter to indicate how many there are
    self.outAtomOrigAtom = {}  # maps to original atom numbers from mol2
    self.outAtomInputConf = {}
    self.outAtomConfNum = {}
    self.confNumAtomList = {}
    for confNum in self.confNums:
      self.confNumAtomList[confNum] = []
      for atomNum in self.confAtoms[confNum]:
        self.outAtoms += 1
        globalAtomNum = self.outAtoms
        self.outAtomOrigAtom[globalAtomNum] = atomNum
        self.outAtomInputConf[globalAtomNum] = self.confInput[confNum][0]
        self.outAtomConfNum[globalAtomNum] = confNum
        self.confNumAtomList[confNum].append(globalAtomNum)

  def _identifyClashSetnums(self, clashDecider, mol2data):
    '''for each set decide if it is broken/clashed
    and add it to the self.brokenConfs list if it is. clashDecider is a
    clash.Clash object that figures out what a clash is. mol2data is the
    mol2.Mol2 object that has atom type information and bondedTo method.'''

    """ 
    Benjamin Tingle 2/13/2021
    Optimized clash sets algorithm for our current one-rule H-H setup
    Only calculates distances between hydrogen atom pairs that match the bond rule
    """
    self.brokenSets = []

    nconformations = len(mol2data.atomXyz)

    #brokenBonds = clashDecider.decideBondRules(mol2data, mol2data.atomXyz[0])
    #if brokenBonds:
    #  print("Bond rule was broken! All conformations are clashed!")
    #  self.brokenSets.extend([aSet for aSet in range(nconformations)])
    #  return

    """for aSet in range(nconformations):
      if clashDecider.decideDistanceRules(mol2data, mol2data.atomXyz[aSet]):
        self.brokenSets.append(aSet)"""

  def _initClusters(self, clusters):
    '''initializes or reinitializes the clusters of conformations'''
    self.clusters = {}
    self.setNameRemap = {}  # maps old sets to new names
    self.setNameOutOrder = []
    self.setNameFirst = {}
    self.setNameLast = {}
    curSetName = 1
    for clusterIndex, cluster in enumerate(clusters):  # save each cluster
      self.clusters[clusterIndex] = tuple(cluster)
      self.setNameFirst[clusterIndex] = curSetName
      for setName in cluster:
        # B.T:
        # it seems this will happen sometimes
        # this is a not good thing to happen so we avoid it now with this
        if self.setNameRemap.get(setName): 
          continue
        self.setNameRemap[setName] = curSetName  # map from old to new
        self.setNameOutOrder.append(setName)
        curSetName += 1  # advance counter
      self.setNameLast[clusterIndex] = curSetName - 1  # doing inclusive

  def _findAdditionalMatchSpheres(self, numSpheres=5, cutoff=2.5):
    '''for each cluster, find a couple matching spheres for distant atoms
    that are relatively localized in space.
    data ends up in dict self.clusterSpheres.
    numSpheres is the max# of spheres to add for each cluster. will not always
     find as many as requested.
    cutoff is used as a cutoff to decide
     whether or not to add a sphere for that atom, mean pairwise dist?'''
    #atomDists = self.mol2data.distFromAtoms(self.rigidComponent)  # useful
    possibleAtoms = set(self.heavyAtomNums)  # only heavy can be matching
    possibleAtoms.difference_update(self.rigidComponent)  # no need to repeat
    #possAtomDist = []  # useful for sorting by distance
    #for possibleAtom in possibleAtoms:
      #possAtomDist.append((possibleAtom, atomDists[possibleAtom]))
    #possAtomDist.sort(key=operator.itemgetter(1), reverse=True)
    #use possAtomDist for each cluster now to find the best candidates
    self.clusterSpheres = {}  # indexed by clusterIndex just like self.clusters
    for clusterIndex in list(self.clusters.keys()):
      cluster = self.clusters[clusterIndex]  # cluster is a tuple of confs
      #print "cluster", cluster #debugging
      self.clusterSpheres[clusterIndex] = []
      for possibleAtom in possibleAtoms:
        xyzPositions = self.mol2data.getXyzManyConfs(cluster, possibleAtom)
        okayToAdd = False
        if 1 == len(xyzPositions):  # singleton cluster, definitely okay
          okayToAdd = True
        else:
          longDist, meanDist = geometry.longestAndMeanDist(xyzPositions)
          if meanDist <= cutoff:  # passes cutoff
            okayToAdd = True
        if okayToAdd:  # either singleton or passes cutoff
          avgPoint = geometry.getAverage(xyzPositions)
          self.clusterSpheres[clusterIndex].append((possibleAtom, avgPoint))
          #print possibleAtom  # debugging
          if len(self.clusterSpheres[clusterIndex]) == numSpheres:  # done
            break  # out of for loop, no need to go on
      #print self.clusterSpheres[clusterIndex] #debugging

  def _makeClouds(self):
    '''highest level of hierachical ligand sampling, breaks the input
    sets into a few clouds representing gross levels of similar conformations'''
    #atomDists = self.mol2data.distFromAtoms(self.rigidComponent)
    #needs switched to divisive bisecting k-means clustering to be fast.
    clusters = self.mol2data.divisiveClustering()
    #printClusterHelper(clusters)  # debug cluster assignments
    self._initClusters(clusters)
    #now that we have clusters, want to find additional matching spheres
    #(with colors even though coloring is bad)
    #data ends up in dict self.clusterSpheres
    self._findAdditionalMatchSpheres()

  def _colorWriter(self, outFile, mol2data):
    '''writes the color table if it was changed from the default'''
    if mol2data.colorConverter.colorInts != \
        mol2data.colorConverter.colorIntsDefault:  # if not default
      colors = list(mol2data.colorConverter.colorInts.items())
      colors.sort(key=operator.itemgetter(1))
      for colorName, colorKey in colors:
        outFile.write('T %2d %8s\n' % (colorKey, colorName))

  def _allButSetWriter(
      self, outFile, mol2data, solvdata, setsTotal, clustersTotal=0):
    '''writes the M A B X R and C lines'''
    #now the molecule section, facts about the whole molecule, 5 lines
    outFile.write(
        'M %16s %9s %3d %3d %6d %6d %6d %6d %6d %6d\n' % (
            mol2data.name[-16:], mol2data.protName[-9:],
            len(mol2data.atomNum), len(mol2data.bondStart),
            self.outAtoms, self.numConfs, setsTotal,
            self.heavyRigidCount, 5, clustersTotal))
    #second molecule line, solvation and charge data
    outFile.write(
        'M %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n' % (
            solvdata.totalCharge, solvdata.totalPolarSolv,
            solvdata.totalApolarSolv, solvdata.totalSolv,
            solvdata.totalSurface))
    #smiles and long version of name
    outFile.write('M %-76s\n' % (mol2data.smiles[-76:]))
    outFile.write('M %-76s\n' % (mol2data.longname[-76:]))
    #best dud energy, computed and put in later. idea is to store the best
    #energy that can be found using the old DOCK/db methods and make sure
    #we aren't totally missing the ball.
    outFile.write('M %+10.4f\n' % 999.999)
    #atom line, 1 per atom
    for atomNum in range(len(mol2data.atomNum)):
      outFile.write(
          'A %3d %-4s %-5s %2d %2d %+9.4f %+10.3f %+10.3f %+10.3f %9.3f\n' % (
              mol2data.atomNum[atomNum], mol2data.atomName[atomNum],
              mol2data.atomType[atomNum],
              mol2data.dockNum[atomNum], mol2data.colorNum[atomNum],
              solvdata.charge[atomNum], solvdata.polarSolv[atomNum],
              solvdata.apolarSolv[atomNum], solvdata.solv[atomNum],
              solvdata.surface[atomNum]))
    #now all the bonds.
    for bondNum in range(len(mol2data.bondStart)):
      outFile.write(
          'B %3d %3d %3d %-2s\n' % (
              mol2data.bondNum[bondNum], mol2data.bondStart[bondNum],
              mol2data.bondEnd[bondNum], mol2data.bondType[bondNum]))
    #now all the coordinates. this section is complex to output since not
    # all atoms*input coordinates are output.
    for xyzNum in range(len(self.outAtomXYZ)):
      #xyzNum += 1  # 1-index nonsense
      atomNum = self.outAtomOrigAtom[xyzNum]
      #inputConfNum = self.outAtomInputConf[xyzNum]
      confNum = self.outAtomConfNum[xyzNum]
      xyz = self.outAtomXYZ[xyzNum]
      #xyz = self.mol2data.atomXyz[inputConfNum][atomNum]
      #atomnum needs incremented by 1 to make it match up with the input atom#
      outFile.write(
          'X %9d %3d %6d %+9.4f %+9.4f %+9.4f\n' %
          (xyzNum+1, atomNum+1, confNum, xyz[0], xyz[1], xyz[2]))
      #amazingly these coordinates are not converted to integers.
    #rigid xyzs, or really just the ligand xyzs to be used for matching
    self.rigidNumSeen = 0
    for rigidNum in self.heavyRigidAtomNums:
      self.rigidNumSeen += 1
      atomColor = mol2data.colorNum[rigidNum]
      xyz = self.mol2data.atomXyz[0][rigidNum]
      outFile.write(
          'R %6d %2d %+9.4f %+9.4f %+9.4f\n' %
          (self.rigidNumSeen, atomColor, xyz[0], xyz[1], xyz[2]))
    #conformations...
    for confNum in range(self.numConfs):
      coordStart = self.confNumAtomList[confNum][0] + 1 # conf start/end should be 1-indexed in the output db2
      coordEnd = self.confNumAtomList[confNum][1] + 1
      outFile.write('C %6d %9d %9d\n' % (confNum+1, coordStart, coordEnd))

  def _setWriter(self, outFile, mol2data, solvdata):
    '''writes the S lines. no more limit here.'''
    #set conf list S
    if self.clusters is None:  # if clusters weren't made
      curSets = list(self.setToConfs.keys())  # this order is fine
      curSets.sort()
    else:
      curSets = self.setNameOutOrder
    for outSetNum, curSet in enumerate(curSets):  # all sets
      if self.clusters is None:  # if clusters weren't made
        outSetNum += 1  # 1 index since it is fortran
      else:
        outSetNum = self.setNameRemap[curSet]
      curConfs = [c+1 for c in self.setToConfs[curSet]]
      totalConfs = len(curConfs)
      if 0 == totalConfs:  # means there are no children, this shouldn't happen
        print("set", curSet, "has no conformations in it.", curConfs)
        sys.exit(1)
      else:
        totalLines = int(math.ceil(totalConfs / self.coSePerLine))
        lastLineLen = totalConfs % int(self.coSePerLine)
        if 0 == lastLineLen:
          lastLineLen += int(self.coSePerLine)  # correct count when 0
        #the first line that says how many more are coming and has data
        inInput = 0  # mix-n-match
        confEnergy_1 = 999999.999
        confEnergy_2 = 999999.999
        outHydro = 3  # mix-n-match
        #this makes the confEnergy a mmff internal energy, ignoring hydroxyls
        #that have been rotated for now.
        #jklyu, 20200511, replace mmff internal energy with totalStrain from UCSF strain energy
        #confEnergy = mol2data.inputEnergy[curSet] - min(mol2data.inputEnergy)
        confEnergy_1 = mol2data.inputTotalStrain[curSet]
        confEnergy_2 = mol2data.inputMaxStrain[curSet]
        #confEnergy = mol2data.inputEnergy[curSet]
        outHydro = mol2data.inputHydrogens[curSet]
        brokenSet = 0  # not broken
        if curSet in self.brokenSets:
          brokenSet = 1  # broken
        outFile.write(
        #    'S %6d %6d %3d %1d %1d %+11.3f\n' % (
        #        outSetNum, totalLines, totalConfs, brokenSet, outHydro,
        #        confEnergy))
            'S %6d %6d %3d %1d %1d %+11.3f %+11.3f\n' % (
                outSetNum, totalLines, totalConfs, brokenSet, outHydro,
                confEnergy_1, confEnergy_2))
        fullLineFormat = 'S %6d %6d %1d'
        for count in range(int(self.coSePerLine)):
          fullLineFormat += ' %6d'
        fullLineFormat += '\n'
        for lineNum in range(totalLines - 1):  # each full line
          outData = [outSetNum, lineNum + 1, self.coSePerLine]
          for count in range(int(self.coSePerLine)):
            outData.append(
                curConfs[lineNum * int(self.coSePerLine) + count])
          outFile.write(fullLineFormat % tuple(outData))
        #now write last line separately and carefully
        partLineFormat = 'S %6d %6d %1d'
        outData = [outSetNum, totalLines, lastLineLen]
        for count in range(lastLineLen):
          partLineFormat += ' %6d'
          outData.append(
              curConfs[(totalLines - 1) * int(self.coSePerLine) + count])
        partLineFormat += '\n'
        outFile.write(partLineFormat % tuple(outData))

  def _cloudWriter(self, outFile, mol2data):
    '''write the cloud data'''
    self.cloudNumSeen = 0
    for clusterId in list(self.clusters.keys()):
      outClusId = clusterId + 1
      countSph = len(self.clusterSpheres[clusterId])
      #next line gets around a bug produced when countSph is 0
      maxSphCount = max(self.cloudNumSeen + countSph, self.cloudNumSeen + 1)
      outFile.write(
          'D %6d %6d %6d %3d %3d %3d\n' % (
              outClusId, self.setNameFirst[clusterId],
              self.setNameLast[clusterId], countSph, self.cloudNumSeen + 1,
              maxSphCount))
      for matchAtom, matchXyz in self.clusterSpheres[clusterId]:
        self.cloudNumSeen += 1  # advance counter
        atomColor = mol2data.colorNum[matchAtom]
        outFile.write(
            'D %6d %2d %+9.4f %+9.4f %+9.4f\n' % (
                self.cloudNumSeen, atomColor,
                matchXyz[0], matchXyz[1], matchXyz[2]))

  def write(
      self, db2gzFileName, verbose=False, timeit=False,
      #limitset=9999999, writeMode='w'):
      limitset=9999999, writeMode='wt'):
    '''writes to the new db2 file format. already gzipped.
    writeMode allows append instead of write(over)'''
    try:  # to open the file
      #outFile = gzip.GzipFile(db2gzFileName, writeMode)
      #outFile = gzip.open(db2gzFileName, writeMode)
      outFile = io.StringIO()
      try:
        mol2data = self.mol2data
      except AttributeError:
        print('mol2data missing when output stage encountered.(3)')
        sys.exit(1)
      try:
        solvdata = self.solvdata
      except AttributeError:
        print('solvdata missing when output stage encountered.(4)')
        sys.exit(1)
      #check if default colors changed, write if they have.
      self._colorWriter(outFile, mol2data)
      setsTotal = len(self.setToConfs.keys())
      self._allButSetWriter(
          outFile, mol2data, solvdata,
          setsTotal, len(self.clusters or []))
      self._setWriter(outFile, mol2data, solvdata)
      if self.clusters is not None:  # if makeclouds was run
        self._cloudWriter(outFile, mol2data)  # this sucks, have to only
           #write clouds for sets that were written. need to rething huge hack
      outFile.write('E\n')  # write the E line here
      outFile.seek(0, 0)
      outgz = gzip.open(db2gzFileName, writeMode)
      outgz.write(outFile.read())
    except IOError:
      print("error opening output file", db2gzFileName)
      sys.exit(1)
    if verbose:
      print(db2gzFileName + " file written out")
    
  def writeFile(
    self, fileHandle, verbose=False, timeit=False,
    limitset=9999999, writeMode='wt'):
    try:
      mol2data = self.mol2data
    except AttributeError:
      print('mol2data missing when output stage encountered.(3)')
      sys.exit(1)
    try:
      solvdata = self.solvdata
    except AttributeError:
      print('solvdata missing when output stage encountered.(4)')
      sys.exit(1)
    outFile=fileHandle
    #check if default colors changed, write if they have.
    self._colorWriter(outFile, mol2data)
    self._allButSetWriter(
        outFile, mol2data, solvdata,
        len(list(self.setToConfs.keys())), 0 if not self.clusters else len(self.clusters))
    self._setWriter(outFile, mol2data, solvdata)
    if self.clusters is not None:  # if makeclouds was run
      self._cloudWriter(outFile, mol2data)  # this sucks, have to only
          #write clouds for sets that were written. need to rething huge hack
    outFile.write('E\n')  # write the E line here

  def writeMol2(
      self, mol2fileName, verbose=False, timeit=False, separateClusters=True):
    '''writes multi-mol2 files instead of db2 files. useful for debugging
    the clustering (or other procedures).     each cluster can be written separately
    and will be given a prefix of cluster.00001. etc'''
    if self.clusters is None:
      separateClusters = False  # don't write non-existent clusters
    if separateClusters:
      currentCluster = list(self.clusters.keys())[0] + 1
      #currentPrefix = "cluster." + string.zfill(currentCluster, 5) + "."
      currentPrefix = "cluster." + currentCluster.zfill(5) + "."
      currentName = currentPrefix + mol2fileName
    else:
      currentName = mol2fileName
    try:  # to open the file
      outFile = open(currentName, 'w')
      try:
        mol2data = self.mol2data
      except AttributeError:
        print('mol2data missing when output stage encountered.(1)')
        sys.exit(1)
      try:
        solvdata = self.solvdata
      except AttributeError:
        print('solvdata missing when output stage encountered.(2)')
        sys.exit(1)
      if self.clusters is not None:  # if makeclouds was run
        outFile.close()  # close the open and empty file. stupid stupid hack.
        for clusterId in list(self.clusters.keys()):
          currentCluster = clusterId + 1
          #currentPrefix = "cluster." + string.zfill(currentCluster, 5) + "."
          currentPrefix = "cluster." + currentCluster.zfill(5) + "."
          currentName = currentPrefix + mol2fileName
          outFile = open(currentName, 'w')
          outNums = []
          for confNumber in range(
              self.setNameFirst[clusterId], self.setNameLast[clusterId] + 1):
            outNum = self.setNameOutOrder[confNumber - 1]  # hate 1-indexing
            outNums.append(outNum)
          self.mol2data.writeMol2File(outFile, outNums)
          if verbose:
            print(currentName + " file written out")
          outFile.close()
      else:
        self.mol2data.writeMol2File(outFile)  # just write them all
      outFile.close()
    except IOError:
      print("error opening output file", currentName)
      sys.exit(1)
    if verbose:
      print(currentName + " file written out")

# fix for py3-3.7
#if -1 != string.find(sys.argv[0], "hierarchy.py"):
if -1 != sys.argv[0].find("hierarchy.py"):
  #nothing to do if called from commandline
  pass
