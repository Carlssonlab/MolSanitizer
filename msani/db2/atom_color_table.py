#!/usr/bin/env python

#Ryan G. Coleman, Brian K. Shoichet Lab
#implements atom color table, 2 parts, 1 part is color to int mapping
# other part is rules for mapping sybyl atoms to colors
# defaults are encoded and writable so users can edit them
# this code is slated for deprecation unless somebody can show that colors help

#import string
import sys
import operator

class ColorTable(object):
  '''implements the color table and all the things it needs'''
  defaultColorDefault = 'neutral'  # the default color
  #note that if you change the 7 classes here you have to change INDOCK
  colorIntsDefault = {'positive':                       1,
                      'negative':                       2,
                      'acceptor':                       3,
                      'donor':                          4,
                      'ester_o':                        5,
                      'amide_o':                        6,
                      'neutral':                        7}
  #rules are last match counts. so even though all Ns are positive, the later
  # rule N.ar matches acceptor so N.ar are acceptor
  #first rules ar  beginning of sybyl atom text -> type
  #later rules are Atom NotBondedTo Atom -> type like
  #                O. -1 N.2 -> negative
  #other rules are Atom BondsAwayFrom Atom -> type like
  #                C.2 1 N. -> positive or
  #                0.2 2 N.3 -> amide_o
  #again rules are read in order and the last matching rule is the one used
  rulesTableDefault = [('N.4', 'positive'),
                       ('O.co2', 'negative'),
                       ('O.2', 'acceptor'),
                       ('O.3', 'acceptor'),
                       ('S.2', 'acceptor'),
                       ('N.ar', 'acceptor'),
                       #('P.3', 1, 'O.co2', 'negative'), #Turned off, as this was chemically incorrect
                       ('O.3', 1, 'P.3', 'negative'), #Added to make O.3 bonded to P.3 negative. Have to fix the ester side later.
                       ('O.3', 1, 'S.o2', 'negative'), #Added to make O.3 bonded to S.o2 negative. Have to fix the ester side later.
                       ('N.am', 1, 'H', 'donor'),
                       ('N.pl3', 1, 'H', 'donor'),
                       #('O.3', 1, 'H', 'donor'), #Duplicated
                       ('N.2', -1, 'H', 'acceptor'),
                       ('N.2', -1, 'C.3', 'acceptor'),
                       ('N.2', 1, 'H', 'donor'),
                       ('N.ar', -1, 'H', 'acceptor'), #These N.ar seems unnecessary but were kept for compatibility with CORINA
                       ('N.ar', -1, 'C.3', 'acceptor'),
                       ('N.ar', 1, 'H', 'donor'),
                       ('O.3', 1, 'H', 'donor'),
                       ('O.2', 2, 'O.3', 'ester_o'),
                       ('O.2', 2, 'N.pl3', 'amide_o'),
                       ('O.2', 2, 'N.am', 'amide_o'),
                       ('O.2', 2, 'N.3', 'amide_o')]

  def __init__(self, parameterFileName=None):
    '''constructs from defaults or reads from file'''
    if parameterFileName is not None:
      parameterFile = open(parameterFileName, 'r')
      phase = 0  # 0 = default, 1 = color/ints, 2 = rules
      self.colorInts = {}
      self.rulesTable = []
      try:
        for line in parameterFile:
          if 0 == phase:
            #self.defaultColor = string.split(line)[0]
            self.defaultColor = line.split()[0]
            phase = 1  # only one line to define default color
          elif 1 == phase:
            #tokens = string.split(line)
            tokens = line.split()
            if 1 == len(tokens) and tokens[0] == 'rules':
              phase = 2  # rules table is next
            else:
              self.colorInts[tokens[0]] = int(tokens[1])
          elif 2 == phase:
            #tokens = string.split(line)
            tokens = line.split()
            if 2 == len(tokens):  # normal rule
              self.rulesTable.append((tokens[0], tokens[1]))
            elif 4 == len(tokens):  # rule about bonds
              self.rulesTable.append(
                  (tokens[0], int(tokens[1]), tokens[2], tokens[3]))
      except StopIteration:
        pass  # EOF
    else:  # no parameter file, use defaults
      self.defaultColor = self.defaultColorDefault
      self.colorInts = self.colorIntsDefault
      self.rulesTable = self.rulesTableDefault

  def printParameters(self):
    '''prints to standard out the parameters used in a readable format'''
    print(self.defaultColor)
    thisColorInts = list(self.colorInts.items())
    thisColorInts.sort(key=operator.itemgetter(1))
    for aColor, anInt in thisColorInts:
      print(aColor, anInt)
    print("rules")
    for rule in self.rulesTable:
      for rulePart in rule:
        print(rulePart, end=' ')
      print("")  # force newline
    #that's all. rules are printed

  def convertMol2color(self, mol2data, atomNum):
    '''uses the rules in the color table to produce the integer of the color'''
    actualNum = atomNum - 1
    actualName = mol2data.atomType[actualNum]
    lastColorFound = self.defaultColor
    #go through each rule in order. apply it if it matches.
    for rule in self.rulesTable:
      if 2 == len(rule):  # easy rule
        #if 0 == string.find(actualName, rule[0]):  # means it matched
        # fix for py3-3.7
        if 0 == actualName.find(rule[0]):  # means it matched
          lastColorFound = rule[1]
      elif 4 == len(rule):  # rule about bonds
        if rule[1] == -1:  # means not bonded to
          # fix for py3-3.7
          if 0 == actualName.find(rule[0]):  # matched the first part
          #if 0 == string.find(actualName, rule[0]):  # matched the first part
            if not mol2data.bondedTo(atomNum, rule[2]):
              lastColorFound = rule[3]
        else:  # has a positive number like 1 or 2
          # fix for py3-3.7
          #if 0 == string.find(actualName, rule[0]):  # matched the first part
          if 0 == actualName.find(rule[0]):  # matched the first part
            if mol2data.bondedTo(atomNum, rule[2], rule[1]):
              lastColorFound = rule[3]
    
    # Special case for ester O.3 bonded to P.3 and C.3 or C.2
    if actualName == 'O.3' and (mol2data.bondedTo(atomNum, 'P.3') or mol2data.bondedTo(atomNum, 'S.o2')) and\
      (mol2data.bondedTo(atomNum, 'C.3') or mol2data.bondedTo(atomNum, 'C.2')):  
      lastColorFound = 'acceptor'
    return self.colorInts[lastColorFound]

# fix for py3-3.7
#if -1 != string.find(sys.argv[0], "atom_color_table.py"):
if -1 != sys.argv[0].find("atom_color_table.py"):
  #if program is called from the command line, assume user wants a copy of the
  #default parameter file written to standard out. this is the only command use.
  #usually this will be imported and run from somewhere else.
  if len(sys.argv) > 1:
    ColorTable(sys.argv[1]).printParameters()
  else:
    ColorTable().printParameters()
