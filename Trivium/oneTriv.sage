from polybori import * 
from utilities import *

#################################################################################
# Polynomial based implementation of Trivum. Needed for the "golden instance"
def triviumGolden(preRounds,R, ivVarNum, **kwds):
  """
  returns the generating quadratic equations and the output of one trivium instance
  
  INPUT:
    R -- polynomial ring. x00..x79: key-variables, i00..i??: ivVariables (at most 80), y??: other
    I -- set of indeces from 0 to 79. for each i in I the function sets the i-th iv bit to 1 (default: [])
    K -- 80 bit key (default: random)
    Nr -- Number of initializing rounds 
    N_o -- number of produced output (must be < 67)
  """
  def update_state(curRound=0):#updates the symbolic state
    #setting new variables
    F.append(a[curRound] + IS[0])
    F.append(b[curRound] + IS[93])
    F.append(c[curRound] + IS[177])
    IS[0] = a[curRound]
    IS[93] = b[curRound]
    IS[177] = c[curRound]
		
    #updating
    t1 = IS[65] + IS[92] + IS[90]*IS[91] + IS[170]
    t2 = IS[161] + IS[176] + IS[174]*IS[175] + IS[263]
    t3 = IS[242] + IS[287] + IS[285]*IS[286] + IS[68]
		
		# shift all registers by 1	
    IS[1:93] = IS[0:92]
    IS[0] = t3
    IS[94:177] = IS[93:176]
    IS[93] = t1
    IS[178:288] = IS[177:287]
    IS[177] = t2

	
  #setting up the BooleanPolynomialRing
  allVars=R.gens()
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  
  iv=[0 for _ in xrange(80)]
  for i in xrange(ivVarNum):
    ivName = "iD_%02d"%i
    iv[i] = R(ivName)
  
  global a
  a = sorted([ v for v in allVars if str(v)[0] == "a"], key=lambda k: int(str(k)[3:]))
  global b
  b = sorted([ v for v in allVars if str(v)[0] == "b"], key=lambda k: int(str(k)[3:]))
  global c
  c = sorted([ v for v in allVars if str(v)[0] == "c"], key=lambda k: int(str(k)[3:]))

  #setting up the state of trivium:
  global F
  F=[]
  global IS
  IS=[R(0) for _ in xrange(288)]
  IS[:80]=x
  IS[93:173]=iv
  IS[285:288] = [R(1),R(1),R(1)]
  #clocking the state for phase1 and phase2
  for curRound in xrange(preRounds):
    update_state(curRound) 
  
  #returning system-equations
  return F

#################################################################################
# Polynomial based implementation of Trivum. Needed for the "golden instance"
def triviumWK(preRounds,R, ivVarNum, keyVarNum, **kwds):
  """
  returns the generating quadratic equations and the output of one trivium instance
  
  INPUT:
    R -- polynomial ring. x00..x79: key-variables, i00..i??: ivVariables (at most 80), y??: other
    I -- set of indeces from 0 to 79. for each i in I the function sets the i-th iv bit to 1 (default: [])
    K -- 80 bit key (default: random)
    Nr -- Number of initializing rounds 
    N_o -- number of produced output (must be < 67)
  """
  def update_state(curRound=0):#updates the symbolic state
    #setting new variables
    F.append(a[curRound] + IS[0])
    IS[0] = a[curRound]
    F.append(b[curRound] + IS[93])
    IS[93] = b[curRound]
    F.append(c[curRound] + IS[177])
    IS[177] = c[curRound]
    
    #updating
    t1 = IS[65] + IS[92] + IS[90]*IS[91] + IS[170]
    t2 = IS[161] + IS[176] + IS[174]*IS[175] + IS[263]
    t3 = IS[242] + IS[287] + IS[285]*IS[286] + IS[68]
    
    # shift all registers by 1  
    IS[1:93] = IS[0:92]
    IS[0] = t3
    IS[94:177] = IS[93:176]
    IS[93] = t1
    IS[178:288] = IS[177:287]
    IS[177] = t2

  #setting up the BooleanPolynomialRing
  allVars=R.gens()
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  
  iv=[0 for _ in xrange(80)]
  for i in xrange(ivVarNum):
    ivName = "iD_%02d"%i
    iv[i] = R(ivName)
  
  global a
  a = sorted([ v for v in allVars if str(v)[0] == "a"], key=lambda k: int(str(k)[3:]))
  global b
  b = sorted([ v for v in allVars if str(v)[0] == "b"], key=lambda k: int(str(k)[3:]))
  global c
  c = sorted([ v for v in allVars if str(v)[0] == "c"], key=lambda k: int(str(k)[3:]))

  #setting up the state of trivium:
  global F
  F=[]
  global IS
  IS=[R(0) for _ in xrange(288)]
  IS[:keyVarNum]=x[:keyVarNum]
  print x[:keyVarNum]
  IS[93:173]=iv
  IS[285:288] = [R(1),R(1),R(1)]
  #clocking the state for phase1 and phase2
  for curRound in xrange(preRounds):
    update_state(curRound) 
  
  #returning system-equations
  return F
##################################################################################################
def outputGolden(startPhase = 150, startOutput = -1, outputBits = -1, outputList = None):
  '''Generates a generic system in the range from/to, plus the corresponding defining
     equations. All variables have the form a/b/c/zD_i with D: instance and i: round number.
     Works its way down from startOutput+outBits towards startPhase.
     Terminate if it finds a full 288 bit state of Trivium and returns the corrected start phase.

     * startPhase: starting point of phase 3 (default: 150). Will be ignored if we 
                   find full 288 state bits in a/b/c before
 
     * startOutput: starting round for output equations. Must be >= startPhase

     * outputBits: number of output bits (minimum: 1)

     * return: goodStart: startPhase identified by the algorithm

     * return: allVars: variables in the system

     * return: outSys: system in Trivium equations and output variables
  '''
  assert startPhase >= 120, "Startphase needs to be >= 120, not %d"%startPhase

  if outputList == None: outputList = [ startOutput..startOutput+outputBits-1 ]
  outputBits = len(outputList)
  startOutput = min(outputList)

  for curRound in outputList:
    assert startPhase <= curRound, "from (%d) must be <= to (%d) with %s"%(startPhase, curRound, outputList)
    assert startOutput <= curRound, "outputlist %s must start at %d, not %d"%(outputList, startOutput, curRound)
  assert len(outputList) > 0, "Need at least one output equation"
  assert startPhase >= 150, "Need to start at 150 or higher, not %d"%startPhase
  
  outSys = []; knownVars = set()
  for curR in outputList:
    # output equation
    i = curR + 1   # funny silly quick fix
    oneEq = [ "zD_%04d"%(curR), "cD_%04d"%(i-66), "cD_%04d"%(i-111), "aD_%04d"%(i-66), "aD_%04d"%(i-93), "bD_%04d"%(i-69), "bD_%04d"%(i-84) ]
    outSys.append("+".join(oneEq))
    knownVars.add("zD_%04d"%(curR))
  # transitively add more defining equations
  keepGoing = True; lastEq = 0; goodStart = None; newVars = set()
  while keepGoing: 
    # extract all variables from the new equations
    allTerms = []
    for curEq in outSys[lastEq:]:
      assert curEq.count(" ") == 0, "No blanks allowed in the equations"
      allTerms.extend(curEq.split("+"))
    lastEq = len(outSys)
    allVars = set(flatten([ t.split("*") for t in allTerms]))
    newVars.update(allVars.difference(knownVars))
    newVars.difference_update(knownVars)
    keepGoing = len(newVars) > 0 # still something to do?

    # check if we have filled all 288 positions
    minA = min([ int(a[3:]) for a in newVars if a[0] == "a" ] + [10000])
    minB = min([ int(b[3:]) for b in newVars if b[0] == "b" ] + [10000])
    minC = min([ int(c[3:]) for c in newVars if c[0] == "c" ] + [10000])
    goodStart = min([ minA+93, minB+84, minC+111 ])
    goodA = filter(lambda a: (goodStart > int(a[3:]) >= goodStart-93) and (a[0] == "a"), newVars)
    goodB = filter(lambda b: (goodStart > int(b[3:]) >= goodStart-84) and (b[0] == "b"), newVars)
    goodC = filter(lambda c: (goodStart > int(c[3:]) >= goodStart-111) and (c[0] == "c"), newVars)
    if ((len(goodA) == 93) and (len(goodB) == 84) and (len(goodC) == 111)) or \
        (len(goodA) == len(goodB) == len(goodC) == 0):
      keepGoing = False
      continue
    else: goodStart = None # we have a problem here...

    maxRound = max([ int(v[3:]) for v in newVars ])
    maxVars = filter(lambda w: int(w[3:]) == maxRound, newVars )
    newVars.difference_update(maxVars)
    # insert transitively all equations we need
    for curVar in maxVars:
      if curVar[0] == "z": continue
      if curVar in knownVars: continue
      i = int(curVar.split("_")[1])
      if i < startPhase: continue

      # let's see what we have to do
      if curVar[0] == "a": oneEq = [ "aD_%04d"%(i), "cD_%04d"%(i-66),"cD_%04d"%(i-111), "cD_%04d*cD_%04d"%(i-110, i-109), "aD_%04d"%(i-69), ]
      if curVar[0] == "b": oneEq = [ "bD_%04d"%(i), "aD_%04d"%(i-66),"aD_%04d"%(i-93), "aD_%04d*aD_%04d"%(i-92, i-91), "bD_%04d"%(i-78), ]
      if curVar[0] == "c": oneEq = [ "cD_%04d"%(i), "bD_%04d"%(i-69),"bD_%04d"%(i-84), "bD_%04d*bD_%04d"%(i-83, i-82), "cD_%04d"%(i-87), ]
      outSys.append("+".join(oneEq))
      # we are done with our work here
      knownVars.add(curVar)

  if goodStart == None: goodStart = startPhase
  
  #delete all z equations
  #outSys = [eq for eq in outSys if eq[0] != 'z']
  #print 'phase3: ',outSys
  
  # extract all variables
  fullVars = set(flatten([ t.split("*") for e in outSys for t in e.split("+") ]))
  
  # we are done here       
  return goodStart, fullVars, outSys


#parallel online-phase output
@parallel
def triviumOutVals(instance, infoList):
  print "instance", instance,
  sys.stdout.flush()
  ivVarNum, varDict, pR, allVars = infoList
  allRounds = [ int(a.split("_")[1]) for a in allVars ]
  minRound = min(allRounds); maxRound = max(allRounds) 
  # set iv variables
  instanceNum = int(instance)
  for i in range(ivVarNum): 
    varDict["iD_%02d"%i] = 1 if ((instanceNum % 2) == 1) else 0
    instanceNum //= 2
  state = initTriviumState(pR, varDict)
  for _ in xrange(minRound):
    updateTrivium(state)
  for r in range(minRound, maxRound+1):
    val = updateTrivium(state)
  return [instance, minRound, maxRound, val]

#parallel data collector
@parallel
def triviumData(instance, infoList):
  print "instance", instance,
  sys.stdout.flush()
  ivVarNum, varDict, allVars = infoList
  allRounds = [ int(a.split("_")[1]) for a in allVars ]
  minRound = min(allRounds); maxRound = max(allRounds) 
  dataEq = []
  # set iv variables
  instanceNum = int(instance)
  for i in range(ivVarNum): 
    varDict["i%02d"%i] = 1 if ((instanceNum % 2) == 1) else 0
    instanceNum //= 2
  
  #init the state of trivium
  fullState = [GF(2)(0) for _ in xrange(288)]
  # setup the key
  for i in xrange(80):
    fullState[i] = GF(2)(varDict['x%02d'%i])
  # setup the IV
  for i in range(ivVarNum):
    fullState[i + 93] = GF(2)(varDict['i%02d'%i])
  # write the final 3 bits
  for i in [285..287]:
    fullState[i] = GF(2)(1)

  for _ in xrange(minRound):
    updateTrivium(fullState)
  for r in range(minRound, maxRound+1):
    if r not in allRounds: updateTrivium(fullState)
    else:
      curVars = [a for a in allVars if int(a.split('_')[1]) == r]
      #appending equations for the state vars
      for curVar in curVars:
        if curVar.startswith('a'): dataEq.append(curVar + '+%d'%fullState[0])
        if curVar.startswith('b'): dataEq.append(curVar + '+%d'%fullState[93])
        if curVar.startswith('c'): dataEq.append(curVar + '+%d'%fullState[177])
      # first part
      tOne = fullState[65] + fullState[92]
      tTwo = fullState[161] + fullState[176]
      tThree = fullState[242] + fullState[287]
      # generate output bit
      zBit = tOne + tTwo + tThree

      #appending equations for the output
      for curVar in curVars:
        if curVar.startswith('z'): dataEq.append(curVar + '+%d'%zBit)
      # second/third part
      tOne = tOne + fullState[90]*fullState[91] + fullState[170]
      tTwo = tTwo + fullState[174]*fullState[175] + fullState[263]
      tThree = tThree + fullState[285]*fullState[286] + fullState[68]

      # shift fullState
      for i in range(287,0,-1):
        if i in [0,93,177]: continue
        fullState[i] = fullState[i-1]
      # set the t-variables
      fullState[0] = tThree
      fullState[93] = tOne
      fullState[177] = tTwo
  
  return dataEq

#################################################################################
# Bit oriented implementation of Trivium. Needed to generate the key stream bits
def initTriviumState(pR, varDict):
  '''Initializes Trivium with a given number of key bits and a start set (IV).

     * pR: PolynomialRing with key variables k00..k79 and iv variables i00..i?? (depending)
  '''

  keyVars = [ pR("x%02d"%a) for a in range(80) ]
  ivVars = [a for a in pR.gens() if str(a)[0] == 'i' ]
  ivVars = [ pR("iD_%02d"%a) for a in range(len(ivVars)) ] 

  # 3 registers:
  #     0..92: key at 0..79
  #    93..176: iv at 93..173 (as much as we have)
  #   177..287: (0..0,1,1,1)

  # generate an empty state
  fullState = []
  for i in range(288):
    fullState.append(pR.zero())

  # setup the key
  for i in xrange(len(keyVars)):
    fullState[i] = keyVars[i]

  # setup the IV
  iv=[0 for _ in xrange(80)]
  fullState[93:173] = iv
  for i in xrange(len(ivVars)):
    fullState[i + 93] = ivVars[i]

  # write the final 3 bits
  for i in [285..287]:
    fullState[i] = pR.one()

  # replace all known variables
  for i in range(288):
    fullState[i] = fullState[i].substitute(varDict)
    assert (fullState[i] == 0) or (fullState[i] == 1), "Only works for constants. Error in bit %d"%i

  return fullState



def updateTrivium(state):
  '''Given a state on individual bits, uses the Trivium update function to generate a new state.
     Clocks exactly one round of Trivium.
  
     * state: state (given from outside)
  '''
  # first part
  tOne = state[65] + state[92]
  tTwo = state[161] + state[176]
  tThree = state[242] + state[287]
  # generate output bit
  zBit = tOne + tTwo + tThree
  # second/third part
  tOne = tOne + state[90]*state[91] + state[170]
  tTwo = tTwo + state[174]*state[175] + state[263]
  tThree = tThree + state[285]*state[286] + state[68]

  # shift state
  for i in range(287,0,-1):
    if i in [0,93,177]: continue
    state[i] = state[i-1]
  # set the t-variables
  state[0] = tThree
  state[93] = tOne
  state[177] = tTwo
  
  return zBit # give the world our result (if needed)


#################################################################################
def specialInstance(pR, ivCounter, eqSys, p3fullSys, rounds, outBits, varValues, ivVarNum, noOutput):
  polyRing = parent(varValues.keys()[0])
  # transfer the iv counter to a setting of the iv variables
  workingVal = ivCounter
  ivDict = {}
  allIVs = ["iD_%02d"%k for k in xrange(ivVarNum)]
  for iv in allIVs:
    ivDict[iv] = str(0)
    varValues[polyRing(iv)] = 0
  # assumption: key values are O.K. We only need to set IV values here
  for curCnt in range(ivVarNum):
    curBit = GF(2)(workingVal % 2)
    #print curCnt+floor(ivCounter/(2**ivVarNum))*ivVarNum
    curIV = "iD_%02d"%curCnt
    varValues[polyRing(curIV)] = curBit
    ivDict[curIV] = str(curBit)
    # finally, go to the next bit
    workingVal //= 2
  print "ivDict", prettyDictStr(ivDict)
  
  workingSys = deepcopy(eqSys)
  p3workingSys = deepcopy(p3fullSys)

  if not(noOutput):
    # generate key stream
    oneState = initTriviumState(pR, varValues)
    for _ in xrange(rounds):
      updateTrivium(oneState)
  
    # add the key equations & save them into the system
    for curRound in range(rounds, rounds+outBits):
      curZ = updateTrivium(oneState)
      #p3workingSys = [ a.replace("zD_%04d"%(curRound), str(curZ)) for a in p3workingSys ]
      p3workingSys.append("z%04d_%04d+%s"%(ivCounter, curRound, curZ))
    print
  # end if

  # replace all IV variables by their corresponding value
  for curIV in ivDict.keys(): # slow but working
    workingSys = [ a.replace(curIV, ivDict[curIV]) for a in workingSys ]
  
  #adding p3 to the system
  workingSys.extend(p3workingSys)

  # replace all "D" by ivCounter (4 digits)
  workingSys = [ s.replace("D", "%04d"%ivCounter) for s in workingSys]

  return workingSys

