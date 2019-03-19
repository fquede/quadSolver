#functions for the trivSolver
from sys import *
import os
import time

# import / attach all necessary files within the project
import configTrivSolver
import satLink, shortEquations, mutants, cnfConvert

from utilities import *

from eqSolver import *
from oneTriv import *
from paintMatrix import *

############################################################################
def instances(args):
  #generating instances and writes the resulting system in a file
  #evaluate which instances should be produced
  ivStart,ivEnd = args.ivBlock.split("/")
  ivStart = int(ivStart); ivEnd = int(ivEnd)
  repIVset = range(ivStart,ivEnd)

  #calculate the number of iv Vars needed to generate the instances
  ivVarNum = 0
  while 2^ivVarNum <= ivEnd: ivVarNum += 1

  #create golden phase3
  p3start, p3vars, p3fullSys = outputGolden(150, args.rounds, args.outBits)

  #ring setup
  stateVars = set(flatten([ ["aD_%04d"%i,"bD_%04d"%i, "cD_%04d"%i] for i in xrange(p3start) ]))
  stateVars.update(p3vars)
  stateVars = sorted(stateVars, key=varname2key)  # sort by round number
  stateVars.reverse()
  varNames = [ "zD_%04d"%j for j in range(args.rounds, args.rounds+args.outBits) ] + stateVars \
             + [ "iD_%02d"%i for i in range(ivVarNum) ] + [ "x%02d"%i for i in range(80) ]
  varNames.reverse()

  pR = BooleanPolynomialRing(len(varNames), varNames, order='degneglex')
  print "start=%d/end=%d"%(ivStart, ivEnd)
  print "number of vars at start: ", len(varNames), "  / start of phase 3", p3start

  # iv variables in a handy dictionary
  valueDict = {}
  for i in range(ivVarNum):
    valueDict[pR("iD_%02d"%i)] = 0 # for the time being - will be changed later

  # create "golden instance"
  goldenEqSys = triviumGolden(p3start, pR, ivVarNum)

  # move the system of equations to an array of strings
  goldenEqSys = [ str(a).replace(" ", "") for a in goldenEqSys ]

  # set up solver
  quadSolver = eqSolver()
  
  # list of new variables for every instance
  varsToQuad = set([ "x%02d"%i for i in range(80) ])
  varsToQuad.update([ "z%04d_%04d"%(i,j) for j in range(args.rounds, args.rounds+args.outBits) for i in range(ivStart, ivEnd) ])
  for v in xrange(ivStart, ivEnd):
    varsToQuad.update([ y.replace("D","%04d"%v ) for y in stateVars ])
  varsToQuad = sorted(varsToQuad, key=varname2key)
  quadSolver.addVariables(varsToQuad)   # define variables for the solver

  nTriv = floor(len(repIVset)/5)+1
  for repCnt in repIVset:
    print "repCnt", repCnt,

    # specialize the instance
    curInstance = specialInstance(pR, repCnt, goldenEqSys, p3fullSys, args.rounds, args.outBits, valueDict, ivVarNum, True)  
    quadSolver.addEquations(curInstance)
  
    if ((repCnt % nTriv) == (nTriv-1)) or (repCnt >= max(repIVset)):
      print
      print "rep", repCnt, "time", time.strftime("%Y-%m-%d %H:%M:%S")
      sys.stdout.flush()

      # start merging
      print "merging the instances"
      evalLoop(quadSolver, args)

  curStat = quadSolver.stats()
  
  f = open(args.writeName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f, "# stats", prettyDictStr(curStat)
  print >>f
  print >>f, "# start lin"
  print >>f, quadSolver.getLinSys()
  print >>f
  print >>f, "# start quad"
  print >>f, quadSolver.getQuadSys()
  f.close()

  # we are done - delete solver
  quadSolver.terminate()

def merge(args):
  #merges systems from different files
  varsToQuad, varsPerFile, fileEqsLin, fileEqsQuad, fileIVStart, fileIVEnd = read_files(args.readFile)
  
  # set up solver
  quadSolver = eqSolver()
  # setting up the order of the solver
  varsToQuad = set([ "x%02d"%i for i in range(80) ]).union(varsToQuad)
  varsToQuad = sorted(varsToQuad, key=varname2key)
  quadSolver.addVariables(varsToQuad)   # define variables for the solver
  
  #storing the order in a file
  orderFileName = args.writeName+'.order'
  f = open(orderFileName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f
  print >>f, "# start order"
  for curVar in varsToQuad:
    print >>f, curVar
  f.close()

  # first: add cubes
  if args.insertCubes != None:
    fullCubeCnt = 0
    print "add cubes",
    sys.stdout.flush()
    cubeEqs = []
    for cubeFile in args.insertCubes:
      print cubeFile,
      sys.stdout.flush()
      f = open(cubeFile, "r")
      for curLine in f:
        # make sure we have some valid input
        if (curLine[0] == "#") or (curLine == "") or (curLine.count('#') == 0) or (curLine.count(':') == 0): continue
        fullCubeCnt += 1
        if (fullCubeCnt % 1000) == 0: print "%d (%d)"%(fullCubeCnt, len(cubeEqs)),   ; sys.stdout.flush()
        curLine = trimStr(curLine)
        # get the important bits
        key, equation = curLine.split(":")
        begin, end = key.split("#")
        if (begin != "shortB") and (end != "shortE") and (not(begin) in varsToQuad) or (not(end) in varsToQuad): continue
        # check if /all/ variables are within varsToQuad
        cubeCnt = equation.count("+")
        if (cubeCnt <= args.cubeLen): cubeEqs.append(equation) # we found a cube :-)
      # end for (curLine)
      f.close() # close current cube file
      print # new line
    # end for cube files
    print "/ total number of cube equations: %d; checked %d"%(len(cubeEqs), fullCubeCnt)
    sys.stdout.flush()
    quadSolver.addEquations(cubeEqs)
  # end insert cubes

  # second: add linear equations
  print "add linear equations"
  for i in range(len(args.readFile)):
    nextFile = args.readFile[i]   
    print nextFile, "at", time.strftime("%Y-%m-%d %H:%M:%S"), "read",
    sys.stdout.flush()
    quadSolver.addEquations(fileEqsLin[i])
    print "/ insert at", time.strftime("%Y-%m-%d %H:%M:%S"),
    sys.stdout.flush()
    evalLoop(quadSolver,args)   # one solving step for integration
    print "/ done at", time.strftime("%Y-%m-%d %H:%M:%S")

  # third: add quadratic equations
  print "add quadratic equations"
  for i in range(len(args.readFile)):
    nextFile = args.readFile[i]   
    print nextFile, "at", time.strftime("%Y-%m-%d %H:%M:%S"), "read",
    sys.stdout.flush()
    quadSolver.addEquations(fileEqsQuad[i])
    print "/ insert at", time.strftime("%Y-%m-%d %H:%M:%S"),
    sys.stdout.flush()
    #evalLoop(quadSolver,args)   # one solving step for integration
    print "/ done at", time.strftime("%Y-%m-%d %H:%M:%S")

  #merge the systems
  evalLoop(quadSolver,args)
  print
  print "finished add"

  # output file
  curStat = quadSolver.stats()
  
  f = open(args.writeName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f, "# stats", prettyDictStr(curStat)
  print >>f
  print >>f, "# start lin"
  print >>f, quadSolver.getLinSys()
  print >>f
  print >>f, "# start quad"
  print >>f, quadSolver.getQuadSys()
  f.close()

  # we are done - delete solver
  quadSolver.terminate()

############################################################################
def export(args):
  #gets a system from strengthen and prepares it for further operations
  #read the system
  varsToQuad, eqSys = read_system(args)

  #get a eqSolver running
  quadSolver=eqSolver()
  quadSolver.addVariables(varsToQuad)
  quadSolver.addEquations(eqSys)

  #echelonize the system
  evalLoop(quadSolver,args)

  #getting the quadratic systems out of the Solver and saving it into a file
  quadSys = quadSolver.getQuadSys().replace(" ", "").split("\n")
  curStat = quadSolver.stats()
  
  f = open(args.writeName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f, "# stats", prettyDictStr(curStat)
  print >>f
  print >>f, "# start lin"
  print >>f, quadSolver.getLinSys()
  print >>f
  print >>f, "# start quad"
  print >>f, quadSolver.getQuadSys()
  f.close()

  # double-counting the variables
  maxVars = []
  #counting per poly
  for poly in quadSys:
    allVars = flatten([term.split('*') for term in poly.split('+')])

    cntDict = { a:0 for a in allVars }
    # get the variable(s) that occurs in the most non-linear monomials
    for curMon in poly.split('+'):
      curMon = curMon.strip()
      if curMon.count('*')>0: # quadratic monom
        for curVar in curMon.split('*'): 
          if cntDict.has_key(curVar): cntDict[curVar] += 1
      elif curMon == '1': continue # no variable 
      else: 
        if cntDict.has_key(curMon): cntDict[curMon] += 1 #linear monomial
    maxVal = max(cntDict.values())
    maxVars.append(filter(lambda v: cntDict[v] == maxVal, cntDict.keys()))
  #counting system-wide
  maxVars = flatten(maxVars)
  varList = {a:maxVars.count(a) for a in set(maxVars)}
  #sorting
  varList = sorted(varList.items(),key=lambda x: x[1],reverse=True)
  if args.guess > 0:
    if not args.guessWrite: args.guessWrite = 'guessVars_R'+str(args.rounds)+'_T' + args.ivBlock.replace("/", "_") + '.txt'
    print "writing guessing vars in File: ", args.guessWrite
    f = open(args.guessWrite,'w')
    f.write('R='+str(args.rounds)+' merge'+'\n')
    f.write(str(varList))
    f.close()
  # done with this method

############################################################################  
def data_collection(args):
  #collect the explicit output and the values for the variables that will be guessed
  #getting data we want to know
  
  varList = []
  if args.minCover:
    print "reading vars to guess out of File: ", args.minCover
    f = open(args.minCover,'r')
    control_string = f.readline()
    print control_string
    varList = eval(f.read())
    f.close()
  

  # extract the iv start and end point
  if args.ivBlock != None: 
    ivStart,ivEnd = args.ivBlock.split("/")
    ivStart = int(ivStart); ivEnd = int(ivEnd)
  
  # create a set of all Triviums that should be generated
  trivInstances = range(ivStart, ivEnd)

  #just in case the varlist is unsorted
  varList = sorted(varList,key=lambda x: x[1])
  varCnt = 0
  fullData = ['z%04d_%04d'%(i,j) for i in trivInstances for j in range(args.rounds,args.rounds+args.outBits)]
  #we are not guessing z-vars (just in case the varList contains some outputVars)
  
  guessData = []
  while varCnt < args.guess:
    curTupel = varList.pop()
    if not(curTupel[0].startswith('z')) and curTupel[0].strip() != '1': guessData.append(curTupel[0].strip()); varCnt += 1
  
  outVals = []
  if args.guessTest:
    for var in guessData:
      curVal = GF(2).random_element()
      print var + '+%d'%curVal
      outVals.append(var + '+%d'%curVal)
  else: fullData+=guessData
  
  # sets the key for the actual instance to break
  if args.keyGen==1: KE = [0 if ((a % 3) == 0) or ((a % 7) == 0) else 1 for a in xrange(80)] #key generation
  else: KE = [GF(2).random_element() for _ in xrange(80)]# random key
  print 'key is: ', KE

  # set number of iv variables needed
  ivVarNum = 0
  while 2^ivVarNum <= ivEnd: ivVarNum += 1

  byInstance = {}
  varDict = { "x%02d"%i: KE[i] for i in range(80) }
  for curVar in fullData:
    if curVar.startswith('x'):
      outVals.append(curVar + '+%d'%varDict[curVar])
      continue
    instance = curVar.split("_")[0][1:]
    if not(byInstance.has_key(instance)): byInstance[instance] = [ivVarNum,varDict,set()]
    byInstance[instance][2].add(curVar)
  print "add data", 
  sys.stdout.flush()

  # generate the corresponding data in parallel
  for X,Y in triviumData(byInstance.items()):
    outVals += Y

  f = open(args.writeName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f, "# the key is:", KE
  print >>f
  print >>f, "# start data"
  for outVal in outVals:
    print >>f, outVal
  f.close()


#######################################################################  
def solveSys(args):
  #solves the system (fast)
  #read the system
  varsToQuad, eqSys = read_system(args)
  
  #read data
  f = open(args.readData, 'r')
  dataEq = []
  for curLine in f:
    if len(curLine) == 0: continue
    if curLine[0] == '#': continue
    dataEq.append(curLine.strip())
  f.close()
  
  #get a eqSolver running
  quadSolver=eqSolver()
  quadSolver.addVariables(varsToQuad)
  quadSolver.addEquations(dataEq)
  quadSolver.addEquations(eqSys)

  print 'start solving the system'
  
  if args.evalSparse: evalNum = int(quadSolver.fullSparseEval()); print "eS", evalNum,
  else: 
    if args.evalMaxOut == None: evalNum = int(quadSolver.fullEval())
    else: evalNum = int(quadSolver.fullEval(outMax=args.evalMaxOut))  
    print "e", evalNum,
  sys.stdout.flush()

  #check consistency
  if args.guessTest:
    knownVars = quadSolver.getVars([ "x%02d"%i for i in range(80) ]).split(',')
    print knownVars
    if knownVars[0][0] == '!': print 'test!'
    else: evalLoop(quadSolver,args)

  #echelonize the system
  if not args.guessTest: evalLoop(quadSolver,args)

  # some statistics
  curStats = quadSolver.stats(level=2)
  print "stats after are", prettyDictStr(curStats)
  sys.stdout.flush()

  # extract key variables
  knownVars = quadSolver.getVars([ "x%02d"%i for i in range(80) ]).split(',')
  print len(knownVars), "known", knownVars

##########################################################################
def check_consistency(args):
  #checks the consistency of the whole system
  #read the system
  varsToQuad, eqSys = read_system(args)
  
  #read data
  f = open(args.readData, 'r')
  dataEq = []
  for curLine in f:
    if len(curLine) == 0: continue
    if curLine[0] == '#': 
      if curLine.startswith('# the key is: '):
        KE = eval(curLine[13:])
        print "key = ", KE
        dataEq+=['x%02d+'%i+str(KE[i]) for i in range(80)]
        print "dataEq:", dataEq
      else: continue
    dataEq.append(curLine.strip())
  f.close()

  #get a eqSolver running
  quadSolver=eqSolver()
  quadSolver.addVariables(varsToQuad)
  quadSolver.addEquations(dataEq)
  quadSolver.addEquations(eqSys)

  print 'start solving the system'

  #echelonize the system
  evalLoop(quadSolver,args)

  # some statistics
  curStats = quadSolver.stats(level=2)
  print "stats after are", prettyDictStr(curStats)
  sys.stdout.flush()

  # extract key variables
  knownVars = quadSolver.getVars([ "x%02d"%i for i in range(80) ]).split(',')
  print len(knownVars), "known", knownVars


#################################################################################
#start help functions
def evalLoop(quadSolver, args):
  lastEval = -1
  while lastEval < 2:
    #main functions
    lastEval += 1
    curStats = quadSolver.stats(level=2)
    # we need echelonize now
    if (int(curStats["quadCols"]) == 0) or (int(curStats["quadRows"]) == 0): break # we have a solution
    #curDensity = 1.0
    #curDensity = 0.0
    curDensity = float(curStats["quadNonZeroEntries"])
    curDensity /= float(curStats["quadRows"])*float(curStats["quadCols"])
    if (curDensity < args.density): quadSolver.echelonizeSparse(); print "S", 
    else: quadSolver.echelonizeDense(); print "D",
    # some extra info for the user
    print "(r=%s / c=%s / v=%s / d=%s)"%(curStats["quadRows"], curStats["quadCols"], curStats["quadVarsInUse"], str(curDensity)),
    sys.stdout.flush()

    #check consistency
    if args.guessTest:
      knownVars = quadSolver.getVars([ "x%02d"%i for i in range(80) ]).split(',')
      print knownVars
      if knownVars[0][0] == '!': print 'test!'; break

    # eval 
    if args.evalSparse: evalNum = int(quadSolver.fullSparseEval()); print "eS", evalNum,
    else: 
      if args.evalMaxOut == None: evalNum = int(quadSolver.fullEval())
      else: evalNum = int(quadSolver.fullEval(outMax=args.evalMaxOut))  
      print "e", evalNum,
    sys.stdout.flush()
    if evalNum != 0: lastEval = 0; continue # keep evaluating 

    #check consistency
    if args.guessTest:
      knownVars = quadSolver.getVars([ "x%02d"%i for i in range(80) ]).split(',')
      print knownVars
      if knownVars[0][0] == '!': print 'test!'; break

    # sl 
    if args.sl >= 0 and lastEval < 2:
      if args.sl == 0: 
        slNum = quadSolver.sl()
        if slNum > 2**63: slNum = -(2**64-slNum)
        print "SL", slNum,
      else: 
        slNum = quadSolver.hsl(curDensity,args.density)
        print "hSL", slNum,
      sys.stdout.flush()
      if slNum != 0: continue  # back to the loop

    print       
  ## end of method

#################################################################################
# save data
outToggle = GF(2)(1)
def writeSystem(args, quadSolver):
  # write the system
  global outToggle
  outToggle += GF(2)(1)
  outName = args.writeName + str(outToggle) 
  print "write to file %s at %s"%(outName, time.strftime("%Y-%m-%d %H:%M:%S"))
  f = open(outName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f
  print >>f, "# start ANF"
  f.write(quadSolver.getQuadSys()); f.write("\n")
  f.write(quadSolver.getLinSys()); f.write("\n")
  f.close()


#################################################################################
# sub-routine of strengthen: starts anf<->cnf converstion
def cnfANFconvert(outerSolver, args):
  # generate a tmp-file name
  fileName = "cnfFile_%s-%s"%(shortNumber(floor(time.time() - 1381220000), 5), os.getpid())
  convObj = cnfConvert.cnfConvert(fileName + ".first.cnf") 

  # fully convert this to our new system
  tmpStats = outerSolver.stats(level=1)
  print "c start output with %s polynomials at %s"%(tmpStats["quadRows"], time.strftime("%Y-%m-%d %H:%M:%S"));  sys.stdout.flush()
  convObj.addList(outerSolver.getQuadSys().split("\n"))

  # start simplificator
  print "c start simplificator at %s"%(time.strftime("%Y-%m-%d %H:%M:%S"));  sys.stdout.flush()
  argList = [cnfConvert.cnfSimplifyProgramme, "-no-bve_gates", "-no-cp3_Abva", "-no-cp3_bva_compl", 
             "-enabled_cp3", "-cp3_stats", "-up", "-subsimp", "-all_strength_res=3", "-no-bva",
             "-bve", "-bve_red_lits=1", "-unhide", "-cp3_uhdIters=10", "-bve_BCElim", "-no-dense"]
  secondFileName = "%s.second.cnf"%(fileName); mapFileName = "%s.second.map"%(fileName)
  argList.append("-dimacs=%s"%secondFileName)  # out file
  argList.append("%s.first.cnf"%(fileName))  # include file name
  cnfSimplifier = subprocess.Popen(args=argList, executable=cnfConvert.cnfSimplifyProgramme, \
                     stdout=subprocess.PIPE, shell=False, \
                     stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
  # collect all output
  retCode = cnfSimplifier.wait() # let the process work its magic
  print "------> begin output cnf simplifier at %s <-------"%(time.strftime("%Y-%m-%d %H:%M:%S"))
  for curLine in cnfSimplifier.stdout:
    print curLine,
  print "------> end output cnf simplifier at %s <-------"%(time.strftime("%Y-%m-%d %H:%M:%S"))
  if retCode != 0: print "Error in simplifier"; return []
  del cnfSimplifier  # delete the process structure

  # delete input file
  os.remove(fileName + ".first.cnf")  
  print "c start reading converted file at %s"%(time.strftime("%Y-%m-%d %H:%M:%S"));  sys.stdout.flush()
  for name in [secondFileName, mapFileName]: 
    print "c read %s at %s"%(name, time.strftime("%Y-%m-%d %H:%M:%S"));  sys.stdout.flush()
    if not(os.path.isfile(name)): print "c   could not find %s"%name; sys.stdout.flush(); continue
    f = file(name, "r")
    buckets = dict()
    for line in f:
      if len(line) == 0: continue
      if (line[0] == "c") or (line[0] == "p"): continue
      curList = line.split(" ")
      curClause = frozenset([ int(a) for a in curList if (a != "") and (int(a) != 0) ])
      if len(curClause) > convObj.maxReadLen: continue
      curBucket = frozenset([ abs(a) for a in curClause ])
      if not(buckets.has_key(curBucket)): buckets[curBucket] = set()
      buckets[curBucket].add(curClause)
    f.close() # we are done reading
    # delete the file after we have read it
    os.remove(name)
  # ??? Can we be a bit more memory friendly here?
  # check the number of buckets
  print "have %d buckets at %s"%(len(buckets), time.strftime("%Y-%m-%d %H:%M:%S"))
  transferCnt = 0; degHigh = 0; outList = []
  collectList = convObj.cnf2xorList(buckets)
  print "did transfer %d buckets at %s"%(len(collectList), time.strftime("%Y-%m-%d %H:%M:%S"))

  # filter non-t results (could be done above for "buckets")
  beforeCnt = len(collectList)
  collectList = filter(lambda e: e.find("t") < 0, collectList)
  print "got %d equations from %d at %s"%(len(collectList), beforeCnt, time.strftime("%Y-%m-%d %H:%M:%S"))
  return collectList


#################################################################################
# main solver loop
def strengthen(args):
    # helper functions
    def mutantXL(args):
      # get everything ready
      mutantEqs = quadSolver.getQuadSys().split("\n")
      mutantSearch = mutants.findMutants(args)
      mutantSum = 0
      while len(mutantEqs) != 0:
        mutantEqs = mutantSearch.getMutants(list(mutantEqs))

        # extract all mutant equations from the corresponding quadSolver
        print "Equations found in last step: ", len(mutantEqs) ; sys.stdout.flush()
        if len(mutantEqs) != 0: 
          mutantSum += len(mutantEqs)
          quadSolver.addEquations(mutantEqs) # add intermediate result
        # store the current snap shot
        writeSystem(args, quadSolver)      
        # we are done here
      return mutantSum

######## End of help functions
    print "before reading: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()

    # read system
    varsToQuad, eqSys = read_system(args)
    quadSolver=eqSolver()
    quadSolver.addVariables(varsToQuad)
    quadSolver.addEquations(eqSys)

    print "after reading: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()
    # first: add cubes
    if args.insertCubes != None:
      fullCubeCnt = 0
      print "add cubes",
      sys.stdout.flush()
      cubeEqs = []
      for cubeFile in args.insertCubes:
        print cubeFile,
        sys.stdout.flush()
        f = open(cubeFile, "r")
        for curLine in f:
          # make sure we have some valid input
          if (curLine[0] == "#") or (curLine == "") or (curLine.count('#') == 0) or (curLine.count(':') == 0): continue
          fullCubeCnt += 1
          if (fullCubeCnt % 1000) == 0: print "%d (%d)"%(fullCubeCnt, len(cubeEqs)),   ; sys.stdout.flush()
          curLine = trimStr(curLine)
          # get the important bits
          key, equation = curLine.split(":")
          begin, end = key.split("#")
          if (begin != "shortB") and (end != "shortE") and (not(begin) in varsToQuad) or (not(end) in varsToQuad): continue
          # check if /all/ variables are within varsToQuad
          cubeCnt = equation.count("+")
          if (cubeCnt <= args.cubeLen): cubeEqs.append(equation) # we found a cube :-)
        # end for (curLine)
        f.close() # close current cube file
        print # new line
      # end for cube files
      print "/ total number of cube equations: %d; checked %d"%(len(cubeEqs), fullCubeCnt)
      sys.stdout.flush()
      quadSolver.addEquations(cubeEqs)
    # end insert cubes

    density = args.density
    # initialize internal variables
    lastEval = -1; outToggle = -1; togMain = True
    while lastEval < 3:
      #main functions
      lastEval += 1

      print "start while: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()

      if togMain or not(args.cnfOn):
        # first, we work with mutants
        numMutEqs = mutantXL(args)
        print "(sum: %d)"%numMutEqs,    ; sys.stdout.flush()
        # when no short equations, sat or mutants stop at lastEval = 2
        if numMutEqs == 0: lastEval += 1
        else: lastEval = 0
        print "after mutantXL: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()
      else: 
        # then, we try cnf/anf-convert
        print "cnfANFconvert: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()
        newEqs = cnfANFconvert(quadSolver, args)
        quadSolver.addEquations(newEqs)
        print "after cnfANFconvert: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()
        if len(newEqs) == 0: lastEval += 1
        else: lastEval = 0
      # use the other one
      togMain = not(togMain)

      curStats = quadSolver.stats(level=2)
      # we need echelonize now
      if (int(curStats["quadCols"]) == 0) or (int(curStats["quadRows"]) == 0): break # we have a solution
      #curDensity = 1.0
      #curDensity = 0.0
      curDensity = float(curStats["quadNonZeroEntries"])
      curDensity /= float(curStats["quadRows"])*float(curStats["quadCols"])
      if (curDensity < density): quadSolver.echelonizeSparse(); print "S", 
      else: quadSolver.echelonizeDense(); print "D",
      sys.stdout.flush()

      print "after echelonization: MEM: %.3f GB (%d byte)"%((1.0*getMemory())/(1000.0^3), getMemory());   sys.stdout.flush()

      # some extra info for the user
      print "(r=%s / c=%s / v=%s)"%(curStats["quadRows"], curStats["quadCols"], curStats["quadVarsInUse"]), 
      sys.stdout.flush()
      
      # eval 
      if args.evalSparse: evalNum = int(quadSolver.fullSparseEval()); print "eS", evalNum,
      else: 
        if args.evalMaxOut == None: evalNum = int(quadSolver.fullEval())
        else: evalNum = int(quadSolver.fullEval(outMax=args.evalMaxOut))  
        print "e", evalNum,
      sys.stdout.flush()
      if (evalNum > 5) or ((evalNum != 0) and int(curStats['quadRows']) <= 20): lastEval = 0; continue # keep evaluating (quick fix for mutants)
      #if evalNum != 0: lastEval = 0; continue # keep evaluating 

      # sl 
      if args.sl >= 0 and lastEval < 2:
        if args.sl==0: slNum = quadSolver.sl()
        else: slNum = quadSolver.hsl(curDensity,density)
        if slNum > 2**63: slNum = -(2**64-slNum)
        if args.sl==0: print "SL", slNum,
        else: print "hSL", slNum,
        sys.stdout.flush()
        if slNum != 0: continue  # back to the loop
      #else: lastEval += 1
      writeSystem(args, quadSolver)      
    ## end loop
  ## end of method

##############################################################################
def read_files(files):
  # read all files (if any), get variables
  fileVars = set(); varsPerFile = []
  fileEqsLin = []; fileEqsQuad = []
  fileIVStart = 0; fileIVEnd = 0
  print "read file(s)",
  sys.stdout.flush()
  for nextFile in files:
    curLinEqs = []; curQuadEqs = []
    newVars = set()
    print nextFile,
    sys.stdout.flush()
    f = open(nextFile, "r")
    for curLine in f:
      if len(curLine) == 0: continue
      if curLine[0] == "#": continue
      curEq = trimStr(curLine)
      if curEq.find("*") >= 0: curQuadEqs.append(curEq)
      else: curLinEqs.append(curEq)
      curVars = getVars(curEq)
      newVars.update(curVars)
    f.close()
    varsPerFile.append(newVars)
    fileVars.update(newVars)
    fileEqsLin.append(curLinEqs); fileEqsQuad.append(curQuadEqs)
  fileVars.difference_update([ "x%02d"%i for i in range(80) ])
  fileIVStart = min([int(a.split('_')[0][1:]) for a in fileVars])
  fileIVEnd = max([int(a.split('_')[0][1:]) for a in fileVars])
  print "start instance %04d, end instance %04d"%(fileIVStart, fileIVEnd)
  print "/ found variables", len(fileVars)
  return [fileVars, varsPerFile, fileEqsLin, fileEqsQuad, fileIVStart, fileIVEnd]
  # end read file variables


############################################################################
def read_system(args):
  # get the order from args.readOrder
  f = open(args.readOrder, 'r')
  varsToQuad = []
  for curLine in f:
    if len(curLine) == 0: continue
    if curLine[0] == '#': continue
    varsToQuad.append(curLine.strip())
  f.close()
  
  #read the system out of args.readFile; mutliple Files are possible File to read!!
  quadEqs = []; linEqs=[]
  for curFile in args.readFile:
    f = open(curFile, 'r')
    for curLine in f:
      if len(curLine) == 0: continue
      if curLine[0] == '#': continue
      curEq = trimStr(curLine)
      if curEq.find("*") >= 0: quadEqs.append(curEq)
      else: linEqs.append(curEq)
    f.close()

  return [varsToQuad, linEqs+quadEqs]

