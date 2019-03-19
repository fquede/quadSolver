# frontEnd of the solver for Trivium
# Use via command line (options see below)

#import python useful stuff
import argparse
from sys import *
import os
import time
import thread

# import / attach all necessary files within the project
import configTrivSolver
import satLink, shortEquations, mutants

from utilities import *

from eqSolver import *
from oneTriv import *
from paintMatrix import *


########################################################################
# get external paramters
parser = argparse.ArgumentParser()

# central switches
parser.add_argument('--ivBlock', default="0/256", type=str,
                    help="Start and end for the block of IV variables. Has the form start/end with start/end being either a number between 0 and 2^80. As for range, the end is *not* included in the list. Default: 0/256")
parser.add_argument('--rounds', '-r', default=560, type=int,
                    help="Number of rounds computed before producing output. Default: 560. Full Trivium uses 1152 rounds")
parser.add_argument('--outBits', '-o', default=1, type=int,
                    help="Number of output rounds used in the equations. Default: 1")

# more central switches
parser.add_argument('--nTriv', '-t', default=192, type=int, 
                    help="Number of trivium instances generated per solver iteration")
parser.add_argument('--keyGen', '-k', default=1, type=int, 
                    help="1 for standard key, 0 for rnd")

#offline-phase
parser.add_argument('--noOutput', action="store_true", 
                    help="Suppresses explicit value for output in the system (default: false).")
parser.add_argument('--fullSysName', default=None, type=str, 
                    help="If set, writes the full system to this file. Overwrites another file of the same name that may already exist.")

#online-phase
parser.add_argument('--insertKey', action="store_true", 
                    help="Add values for the key variables x00..x79.")
parser.add_argument('--insertOutput', action="store_true", 
                    help="Add values for the output variables zD_r (D: instance, r: round).")
parser.add_argument('--readFile', '-f', action="append", 
                    help="Uses a given file as auxiliary input. Multiple files are possible.")
parser.add_argument('--readOrder', '-ro', action="append", 
                    help="Uses a given file as order")
parser.add_argument('--ivEnd', '-fe', default=0, type=int, 
                    help="Choose when to start with generating instances")
parser.add_argument('--justFiles', action="store_true", 
                    help="Do not generate additional Trivium instances for solving.")


# parameters for quadSolver
parser.add_argument('--evalMaxOut', default=None, type=int, 
                    help="Maximal number of monomials when using eval")
parser.add_argument('--density', '-d', default=0.0016, type=float, 
                    help="value between 0 and 1. Gives the minimum density for which we use M4RI")
parser.add_argument('--sl', '-s', default=0, type=int, 
                    help="0 for sparse. 1 for half sparse. -1 disables SL")
parser.add_argument('--eval', '-e', default=1, type=int, 
                    help="0 for sparse. 1 for dense")

# sat
parser.add_argument('--satSolver', action="store_true", 
                    help="Enables SAT-solving. Uses CryptoMiniSat as external tool. Default: False")
parser.add_argument('--satKeepFiles', action="store_true", 
                    help="keep intermediate files in case of SAT-solving. Default: False")
parser.add_argument('--satMaxTime', default=60, type=int, 
                    help="Maximal time for one invocation of the SAT-solver. Unit: Seconds. Default: None (no time limit)")
parser.add_argument('--satMaxRestarts', default=None, type=int, 
                    help="Maximal number of restarts for the SAT-solver. Default: None (maximal number of restarts)")
parser.add_argument('--satSplitSize', default=-1, type=int, 
                    help="Maximal size of an equation that we attempt to split.")
parser.add_argument('--satDirectConversion', default=6, type=int, 
                    help="Maximal number of variables for direct conversion between ANF (Xor-Solver) and CNF (SAT-Solver). Default: 6")

# short equations
parser.add_argument('--satShort', action="store_true", 
                    help="Creates short equations using a Genetic Algorithm. Default: False")
parser.add_argument('--satShortMaxStore', default=12, type=int, 
                    help="Maximal size of an equation that we store in the GA.")

parser.add_argument('--directShort', action="store_true", 
                    help="Creates short equations of degree 2. Directly inserts them into the quadSolver. Default: False")
parser.add_argument('--directShortMaxPolys', default=5000, type=int, 
                    help="Maximal number of polynomials considered in direct short equations of degree 2 (default: 5000)")
parser.add_argument('--directShortMaxNewMons', default=1, type=int, 
                    help="Maximal number of new monomials per polynomial (default: 1)")
parser.add_argument('--directShortMaxMem', default=4, type=float, 
                    help="Maximal amount of memory to generate short equations (default: 4, minimum: 1)")
parser.add_argument('--directShortMaxLen', default=6, type=int, 
                    help="Maximal length of polynomials considered for short equations (default: 6)")


# Mutants 
parser.add_argument('--mutantSoftTerms', default=0, type=int, 
                    help="Maximal number of terms that are considered for degree 3 mutants. Values of 0 or less disables soft mutants. Sensible choices are 10..20. Default: 0.")
parser.add_argument('--mutantSoftMaxWork', default=20000, type=int, 
                    help="Maximal number of multiplilcations done for soft degree 3 mutants. Default: 20000.")
parser.add_argument('--mutantSoftMaxCopy', default=20000, type=int, 
                    help="Maximal number of polynomials copied for soft degree 3 mutants. Default: 20000.")
parser.add_argument('--mutantSoftMaxExpand', default=2, type=int, 
                    help="Maximal number of new terms that are introduced by each new polynomial. Default: 2.")

parser.add_argument('--mutantHardTerms', default=0, type=int, 
                    help="Maximal number of terms that are considered for degree 3 mutants. Values of 0 or less disables hard mutants. Sensible choices are 5..20. Default: 0.")
parser.add_argument('--mutantHardMaxCopy', default=5000, type=int, 
                    help="Maximal number of polynomials copied for hard degree 3 mutants. Default: 5000.")

parser.add_argument('--mgaOn', action="store_true", 
                    help="Activates Genetic Algorithms (GA) for Mutants. Requires --mutantTerms > 0. Default: False")
parser.add_argument('--mgaInput', default=2000, type=int, 
                    help="Maximal number of polynomials that are used as input for the GA. Default: 2000.")
parser.add_argument('--mgaRepeat', default=20, type=int, 
                    help="Maximal number of (internal) repetitions of the GA. Default: 20.")
parser.add_argument('--mgaTerms', default=30, type=int, 
                    help="Maximal number of terms for the polynomials in the GA. Default: 30.")
parser.add_argument('--mgaPolys', default=500, type=int, 
                    help="Maximal polynomials considered in the GA. Default: 500.")


#guessing
parser.add_argument('--guess', '-g', default=0, type=int, 
                    help="Guessing up to n variables. 0 for no guessing. Should be less than 40. (default: 0)")
parser.add_argument('--minCover', '-mc', default=None, type=str, 
                    help="sets variables to guess from a file. If None finds variables to guess and writes them to a file (default:None)")
parser.add_argument('--guessRepCnt', '-grc', default=0, type=int,
                    help="0 for guessing at the first solver Iteration. n for guessing at repCnt >= guessrepCnt. There will be just one guessing step! (default: 0)")
parser.add_argument('--guessWrite', '-gw', default=None, type=str,
                    help="writes the variables to guess in a file (see below for default)")
# other parameters
parser.add_argument('--paint', '-a', default=0, type=int, 
                    help="0 for no matrix painting. For options see the comments at the function")
parser.add_argument('--writeSystem', default=-1, type=int, 
                    help="interval in which the system is written (linear and quadratic equations). Default: -1")
parser.add_argument('--insertCubes', action="append", 
                    help="Adds additional cubes equations. Multiple files are possible")

# verbose
parser.add_argument('--verbose', action="append", 
                    help="add verbosity parameters (seperated by comma)")


########################################################################


#evalute external parameters
args = parser.parse_args()
global outBits
outBits = args.outBits
rounds = args.rounds
nTriv = args.nTriv
ckey = args.keyGen
density = args.density
slMode = args.sl
evalMode = args.eval
guessNum = args.guess
guessFlag = guessNum > 0
gRepCnt = args.guessRepCnt
paint = args.paint
haveFiles = args.readFile != None

# check if the paramters are within their range
assert outBits <= 66, "Code only supports 66 outputBits for Trivum."
#if (args.roundsUntil != None) and (args.roundsFrom != None):
#  assert args.roundsUntil <= roundsFrom, "Empty round filter until=%d > from=%d"%(args.roundsUntil, args,roundsFrom)
if haveFiles: 
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file"%curFile
if args.insertCubes != None: 
  for curFile in args.insertCubes:
    assert os.path.isfile(curFile), "%s is not a file"%curFile
assert not(args.mgaOn) or (args.mutantHardTerms > 0) or (args.mutantSoftTerms > 0), "Need mutantTerms > 0 if mgaOn"

# verbosity set
if args.verbose != None: 
  verbStr = ",".join(args.verbose)
  verbositySet = set(verbStr.split(","))
else: verbositySet = set()

# set some internal variables
knownX = dict()
repCntInc = 0
internal = {}
repCnt = 0 

#args.toGuess=range(22,80)
#print args.toGuess,len(args.toGuess)

# needed for correct writing of intermediate results
writeSysCnt = 0

#################################################################################
# main solver loop
def solverLoop(quadSolver):
    def shortSATeq():
      #command for short equations
      print "shortSAT", 
      shortWideEqs.addEquations(quadSolver.getQuadSys())
      newShort = shortWideEqs.createNewEquations(resSize=500) 
      # !!! transfer to SAT-solver
      # filter out quadratic equations, add to quadSolver
      fineVars = getVars(quadSolver.getQuadSys().replace("\n", "+").replace(" ", ""))
      addList = filter(lambda c: max([ a.count("*")+1 for a in c.split("+") ]) <= 2, newShort)
      addList = filter(lambda d: getVars(d).issubset(fineVars), addList)

      # add equations to quadSolver
      quadSolver.addEquations(addList)
      return len(addList), len(newShort)
      
    def mutantXL():
      if (args.mutantHardTerms <= 0) and (args.mutantSoftTerms <= 0): return 0  # nothing to do 
      # get everything ready
      mutantEqs = quadSolver.getQuadSys()
      mutantSum = 0; again = True
      mutRows = -1
      while again:
        again = False
        if args.mutantHardTerms > 0: 
          print "hM", ; sys.stdout.flush()
          mutantSearch.addHardMutants(mutantEqs, maxPolyLen=args.mutantHardTerms, maxCopy=args.mutantHardMaxCopy, \
            verbose=("mutants" in verbositySet))
        if args.mutantSoftTerms > 0: 
          print "sM", ; sys.stdout.flush()
          mutantSearch.addSoftMutants(mutantEqs, maxTerms=args.mutantSoftTerms, maxMuls=args.mutantSoftMaxWork, \
            maxCopy=args.mutantSoftMaxCopy, maxExpand=args.mutantSoftMaxExpand, verbose=("mutants" in verbositySet))

        # extract all mutant equations from the corresponding quadSolver
        mutantEqs = mutantSearch.extractMutants(verbose=("mutants" in verbositySet))
        print len(mutantEqs), ; sys.stdout.flush()
        if len(mutantEqs) != 0: 
          mutantSum+=len(mutantEqs)
          quadSolver.addEquations(mutantEqs) # add intermediate result
          if (args.mgaOn) and (len(mutantEqs) > 50): again = True   # another round
          if not(args.mgaOn): again = True

        if args.mgaOn:
          # also add GA-mutants
          print "gM",
          sys.stdout.flush()
          newEqs = mutantSearch.getInsert(inSys=quadSolver.getQuadSys(), maxCopy=args.mgaInput, repeatNum=args.mgaRepeat, maxPolyLen=args.mgaTerms, \
                     resSize=args.mgaPolys, verbose=("mutants" in verbositySet) )
          print len(newEqs),
          sys.stdout.flush()
          if len(newEqs) != 0: 
            mutantSum+=len(newEqs)
            quadSolver.addEquations(newEqs) # add intermediate result

        # we are done here
      return mutantSum

    def satTry():
      #command for the SAT solver
      print "SAT",
      sys.stdout.flush()
      satSolver.addList(quadSolver.getQuadSys())
      satSolver.addList(quadSolver.getLinSys())
      quadEqs = satSolver.solve(maxTime=args.satMaxTime, maxRestarts=args.satMaxRestarts, keepFiles=args.satKeepFiles)
      print "%d"%(len(quadEqs)),
      sys.stdout.flush()
      if args.satShort: shortWideEqs.addEquations(quadEqs)
      quadSolver.addEquations(quadEqs)
      return len(quadEqs)

######## End of help functions
    lastEval = -1
    solverStep = 0
    doShortCnt = 1
    satLoop = False
    if args.satSolver: satLoop = True
    while lastEval < 3:
      #main functions
      solverStep += 1
      lastEval += 1
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

      # some extra info for the user
      print "(r=%s / c=%s / v=%s)"%(curStats["quadRows"], curStats["quadCols"], curStats["quadVarsInUse"]),
      
      # paint the matrix (if required)
      if paint != 0: 
        namePic = "r" + str(rounds) + "orderRSI_" + "rep" + str(repCnt) + "step" + str(solverStep)
        paintPic(paint, quadSolver.getQuadSys(), namePic)
  
      # write the system (if required)
      if (args.writeSystem != -1):
        nameWrite = "fullSys_r" + str(rounds) + "o" + str(outBits) + "rep" + str(repCnt) + "step" + str(solverStep) + ".anf"
        f = open(nameWrite, "w")
        f.write(quadSolver.getQuadSys()); f.write("\n")
        f.write(quadSolver.getLinSys()); f.write("\n")
        f.close()

      # eval 
      if evalMode == 0: evalNum = int(quadSolver.fullSparseEval()); print "eS", evalNum,
      else: 
        if args.evalMaxOut == None: evalNum = int(quadSolver.fullEval())
        else: evalNum = int(quadSolver.fullEval(outMax=args.evalMaxOut))  
        print "e", evalNum,
      sys.stdout.flush()
      #if (evalNum > 20) or (evalNum != 0 and int(curStats['quadRows']) <= 20): lastEval = 0; continue # keep evaluating (quick fix for mutants)
      if evalNum!=0: lastEval = 0; continue # keep evaluating 

      # sl 
      if slMode >= 0 and lastEval < 2:
        if slMode==0: slNum = quadSolver.sl()
        else: slNum = quadSolver.hsl(curDensity,density)
        if slNum > 2**63: slNum = -(2**64-slNum)
        if slMode==0: print "SL", slNum,
        else: print "hSL", slNum,
        sys.stdout.flush()
        if slNum != 0: continue  # back to the loop
      #else: lastEval += 1

      # do /all/ of these steps
      if args.directShort and lastEval < 3:
        shortList = filter(lambda l: l.count("+") <= args.directShortMaxLen, quadSolver.getQuadSys().split("\n"))
        allowedMons = set([ trimStr(a) for b in quadSolver.getQuadSys().split("\n") for a in b.split("+") if a.count("*") == 1 ])
        print "sd",   ; sys.stdout.flush()
        newEqs = shortMaker.getMore(shortList, allowedMons, args.directShortMaxNewMons, args.directShortMaxPolys)
        print len(newEqs),   ; sys.stdout.flush() 
        quadSolver.addEquations(newEqs)
      if args.satShort and lastEval < 3:
        quadNew,allNew = shortSATeq()
        print "%d/%d"%(quadNew,allNew),   ; sys.stdout.flush()
      if satLoop:
        satLoop = False # only once!
        numSatEqs = satTry()
        print numSatEqs,  ; sys.stdout.flush()
      if ((args.mutantHardTerms != 0) or (args.mutantSoftTerms != 0)) and (lastEval < 3):
        numMutEqs = mutantXL()
        print "(sum: %d)"%numMutEqs,    ; sys.stdout.flush()

      # when no short equations, sat or mutants stop at lastEval = 2
      if not(args.directShort) and not(args.satShort) and not(satLoop) and \
        (args.mutantSoftTerms == 0) and (args.mutantHardTerms == 0): lastEval += 1
      
    print       
  ## end of method
  

#################################################################################
# start of the main programme
#################################################################################
# tell the world what we are doing
startTime=str(time.strftime("%Y-%m-%d %H:%M:%S"))
print "start at", time.strftime("%Y-%m-%d %H:%M:%S")
print "param string", sys.argv, " / all parameters are (line below)"
print prettyDictStr(vars(args))

print
sys.stdout.flush()

#################################################################################

# get the order from args.readOrder
fileIVStart = 0; fileIVEnd = 0
varsToQuad = set()
if args.readOrder:
  for orderFile in args.readOrder:
    f = open(orderFile, 'r')
    for curLine in f:
      if len(curLine.strip()) == 0: continue
      if curLine[0] == '#': continue
      varsToQuad.add(curLine.strip())
    f.close()
  varsToQuad.difference_update([ "x%02d"%i for i in range(80) ] + ['k%03d'%i for i in xrange(600)])

  fileIVStart = int(min([ var.split('_')[0][1:] for var in varsToQuad ]))
  if args.ivEnd != 0: fileIVEnd = args.ivEnd
  else: fileIVEnd = int(max([ var.split('_')[0][1:] for var in varsToQuad ]))
print "start/end of katan instances out of Files", fileIVStart, "/", fileIVEnd

#read the system out of args.readFile
quadEqs = []; linEqs=[]
if args.readFile:
  for dataFile in args.readFile:
    f = open(dataFile, 'r')
    for curLine in f:
      if len(curLine.strip()) == 0: continue
      if curLine[0] == '#': continue
      curEq = trimStr(curLine)
      if curEq.find("*") >= 0: quadEqs.append(curEq)
      else: linEqs.append(curEq)
    f.close()

if args.readFile:
  print 'done reading %d file(s)'%len(args.readFile)
  sys.stdout.flush()

# conversion of all IV variables
# extract the iv start and end point
if args.ivBlock != None: 
  ivStart,ivEnd = args.ivBlock.split("/")
  print ivStart
  ivStart = int(ivStart); ivEnd = int(ivEnd)

# create a set of all Triviums that should be generated
if fileIVEnd!=0: fileIVEnd+=1
knownTrivium = set(range(fileIVStart,fileIVEnd))
if args.justFiles: wantedTrivium = set()
else: wantedTrivium = set(range(ivStart, ivEnd))

repIVset = sorted(wantedTrivium.difference(knownTrivium))
# check how many instances we have in our current system (including files)
curEnd = max(fileIVEnd,ivEnd)
ivVarNum = 0
while 2^ivVarNum <= curEnd: ivVarNum += 1

# sets the key for the actual instance to break
if ckey==1: KE = [0 if ((a % 3) == 0) or ((a % 7) == 0) else 1 for a in xrange(80)] #key generation
else: KE = [GF(2).random_element() for _ in xrange(80)]# random key
print "the key is: ", KE

#create golden phase3
p3start, p3vars, p3fullSys = outputGolden(150, rounds, outBits)

#ring setup
stateVars = set(flatten([ ["aD_%04d"%i,"bD_%04d"%i, "cD_%04d"%i] for i in xrange(p3start) ]))
stateVars.update(p3vars)
stateVars = sorted(stateVars, key=varname2key)  # sort by round number
stateVars.reverse()
varNames = [ "zD_%04d"%j for j in range(rounds, rounds+outBits) ] + stateVars + [ "iD_%02d"%i for i in range(ivVarNum) ] + [ "x%02d"%i for i in range(80) ]
varNames.reverse()

pR = BooleanPolynomialRing(len(varNames), varNames, order='degneglex')
print "start=%d/end=%d"%(ivStart, ivEnd)
print "number of vars at start: ", len(varNames), "  / start of phase 3", p3start

# iv & key variables in a handy dictionary
valueDict = {}
for i in range(80):
  valueDict[pR("x%02d"%i)] = KE[i]
for i in range(ivVarNum):
  valueDict[pR("iD_%02d"%i)] = 0 # for the time being - will be changed later

# create "golden instance"
goldenEqSys = triviumGolden(p3start, pR, ivVarNum)

# move the system of equations to an array of strings
goldenEqSys = [ str(a).replace(" ", "") for a in goldenEqSys ]

#list of key variables
keyVars = [ "x%02d"%i for i in range(80) ]

# set up solver
quadSolver = eqSolver()

# list of new variables for every instance
varsToQuad.update(keyVars)
varsToQuad.update([ "z%04d_%04d"%(i,j) for j in range(rounds, rounds+outBits) for i in range(ivStart, ivEnd) ])
for v in xrange(ivStart, ivEnd):
  varsToQuad.update([ y.replace("D","%04d"%v ) for y in stateVars ])
varsToQuad = sorted(varsToQuad, key=varname2key)
quadSolver.addVariables(varsToQuad)   # define variables for the solver


#storing the order in a file
if args.fullSysName:
  orderFileName = args.fullSysName+'.order'
  f = open(orderFileName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f
  print >>f, "# start order"
  for curVar in varsToQuad:
    print >>f, curVar
  f.close()

#set up the solver for guessing
if guessFlag:
  quadCheater = eqSolver()
  quadCheater.addVariables(varsToQuad)
  if gRepCnt == 0: gRepCnt = ivEnd


varsToQuad = set(varsToQuad)        # for later usage

#adding equations for the known vars in the quadSolver
#quadSolver.addEquations(['x%02d+'%i + str(valueDict[pR('x%02d'%i)]) for i in args.toGuess])

#############################################################################################
# add key values
if args.insertKey:
  keyVals = [ "x%02d+%d"%(i, KE[i]) for i in range(80) ]
  quadSolver.addEquations(keyVals)
  print "added key", keyVals

#############################################################################################
# add output values
if args.insertOutput:
  fullOutput = filter(lambda v: v[0] == "z", varsToQuad)
  byInstance = {}
  varDict = { "x%02d"%i: KE[i] for i in range(80) }
  for curVar in fullOutput:
    instance = curVar.split("_")[0][1:]
    if not(byInstance.has_key(instance)): byInstance[instance] = [ivVarNum,varDict,pR,set()]
    byInstance[instance][3].add(curVar)
  print "add output", 
  sys.stdout.flush()
  # set key variables
  outVals = []
  # generate the corresponding output bit in parallel
  for X,Y in triviumOutVals(byInstance.items()):
    instance, minRound, maxRound, val = Y
    for r in range(minRound, maxRound+1):
      outVals.append('z'+ instance + "_%04d+%d"%(r,val) )
  # add all output values
  quadSolver.addEquations(outVals)
  print "--- added %d values for min/max=%d/%d"%(len(outVals), minRound, maxRound)


#############################################################################################
# add cubes
if args.insertCubes != None:
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
      curLine = trimStr(curLine)
      # get the important bits
      key, equation = curLine.split(":")
      begin, end = key.split("#")
      if (begin != "shortB") and (end != "shortE") and (not(begin) in varsToQuad) or (not(end) in varsToQuad): continue
      # check if /all/ variables are within varsToQuad
      skipCube = False
      cubeVars = getVars(equation)
      if not(cubeVars.issubset(varsToQuad)): continue
      if haveFiles:
        cnt = -1
        for fileVars in varsPerFile: 
          cnt += 1
          # was the cube already in a sub-file?
          if (begin in fileVars) and (end in fileVars): skipCube = True; break
        # end for fileVars  
      if skipCube: continue  # cube already covered in a file
      # finally, we found a cube :-)
      cubeEqs.append(equation)
      sys.stdout.flush()
    # end for (curLine)
    f.close()
  # end for cube files
  print "/ total number of cube equations", len(cubeEqs)
  sys.stdout.flush()
  quadSolver.addEquations(cubeEqs)
# end insert cubes

# enable the sat solver (if applicable)
if args.satSolver:
  satSolver = satLink.satXorWrapper(_splitSize = args.satSplitSize, _cnfMaxLen=args.satDirectConversion)        
if args.satShort:
  shortWideEqs = shortEquations.gaShortWide(_maxVarsReturn = args.satDirectConversion)
  doShortCnt = 0
if args.directShort:
  shortMaker = shortEquations.shortDegTwo()
if (args.mutantSoftTerms > 0) or (args.mutantHardTerms > 0):
  mutantSearch = mutants.findMutants()


#add equations from Files to the quadSolver
if args.readFile:
  quadSolver.addEquations(linEqs)
  quadSolver.addEquations(quadEqs)

## debug
#    print
  print "finished add"
# end read file variables
#############################################################################################

#######################
# delete variables
del varsToQuad, linEqs, quadEqs
#######################



# set the maximal number of rounds to terminate the algorithm
trivCount = nTriv

#if justFiles or repIVset is empty make sure that there is one full iteration
if repIVset == []: repIVset.append('oI')
print repIVset
for repCnt in repIVset:
  print "repCnt", repCnt,

  # specialize the instance
  if repCnt != 'oI':
    curInstance = specialInstance(pR, repCnt, goldenEqSys, p3fullSys, rounds, outBits, valueDict, ivVarNum, args.noOutput)  
    quadSolver.addEquations(curInstance)
  
  # start solving - if we have enough instances (or are in the last round anyway)
  if repCnt == 'oI': repCnt = gRepCnt + 1 #make sure we make a full iteration
  
  if ((repCnt % trivCount) == (trivCount-1)) or (repCnt >= max(repIVset)):
    print
    print "rep", repCnt, "time", time.strftime("%Y-%m-%d %H:%M:%S")
    print "   stats before are", prettyDictStr(quadSolver.stats(level=2))
    sys.stdout.flush()

    #setting trivCount (big step / baby step)
    if (trivCount%nTriv) == 0: trivCount = floor(nTriv/4)
    else: trivCount = nTriv

    # start solving
    #solverLoop(quadSolver)
    if not(guessFlag) or (repCnt < gRepCnt-1): solverLoop(quadSolver)
    
    # check if we need to write the system
    if args.fullSysName != None:
      print "write full system to file %s"%args.fullSysName
      
      f = open(args.fullSysName, "w")
      print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
      print >>f, "# parameter string was:", sys.argv
      print >>f, "# all parameters were:", prettyDictStr(vars(args))
      print >>f, "# stats", prettyDictStr(quadSolver.stats())
      print >>f
      print >>f, "# start lin"
      print >>f, quadSolver.getLinSys()
      print >>f
      print >>f, "# start quad"
      print >>f, quadSolver.getQuadSys()
      f.close()
    
    # some statistics
    curStats = quadSolver.stats(level=2)
    print "stats after are", prettyDictStr(curStats)
    sys.stdout.flush()

    # guess in the last round, after the first full iteration of the solver
    if guessFlag and (repCnt >= gRepCnt-1):
      guessFlag = False  # only guess once
      startGuessTime=str(time.strftime("%Y-%m-%d %H:%M:%S"))
      print 'start guessing at ',startGuessTime


      #getting the quadratic systems out of the Solver
      quadSys = quadSolver.getQuadSys().replace(" ", "").split("\n")
      quadSys.pop() # get rid of the zero
      #getting the variables to guess
      if args.minCover:
        print "reading vars to guess out of File: ", args.minCover
        f = open(args.minCover,'r')
        control_string = f.readline()
        print control_string
        control_string = control_string.split(',')
        assert (control_string[0][2:].strip() == str(rounds) and control_string[1][8:].strip() <= str(gRepCnt)), \
                "variables to guess are not for the current problem"
        varList = eval(f.read())
        f.close()
      else:
        # double-counting the variables
        maxVars = []
        #counting per poly
        for poly in quadSys:
          allVars = flatten([term.split('*') for term in poly.split('+')])

          cntDict = { a:0 for a in allVars}
          # get the variable(s) that occurs in the most non-linear monomials
          for curMon in poly.split('+'):
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
        if not args.guessWrite: args.guessWrite = 'guessVars_R'+str(rounds)+'_T'+str(repCnt) + '.txt'
        print "writing vars in File: ", args.guessWrite
        f = open(args.guessWrite,'w')
        f.write('R='+str(rounds)+',gRepCnt='+str(gRepCnt)+'\n')
        f.write(str(varList))
        f.close()

      if not(args.noOutput):
        #we are in the online phase
        #just in case the varlist is unsorted
        varList = sorted(varList,key=lambda x: x[1])
        varCnt = 0
        varToGuess = []
        #we are not guessing z-vars (just in case the varList contains some outputVars)
        while varCnt < guessNum:
          curTupel = varList.pop()
          if not(curTupel[0].startswith('z')): varToGuess.append(curTupel[0]); varCnt += 1
        print varToGuess
        
        #finding the value of all variables in varList.keys()
        linSys = quadSolver.getLinSys().replace(" ", "").split("\n")
        quadCheater.addEquations(linSys)
        quadCheater.addEquations(quadSys)
        for i in xrange(80):
          quadCheater.addEquations(["x%02d"%i + "+" + str(KE[i])])
    
        #solving the cheated System
        print "quadCheater loop"
        solverLoop(quadCheater)
        
        guessedVars = quadCheater.getVars(varToGuess).split(",")
        #adding the guessed values to the Solver
        print "x are %d from"%len(filter(lambda v: v[0] == "x", guessedVars)), guessedVars 
        quadSolver.addEquations([s.replace("=", "+") for s in guessedVars])     

        print "we are done guessing"
        #kill the quadCheater
        quadCheater.terminate()

        #time used for guessing 
        endGuessTime = time.strftime("%Y-%m-%d %H:%M:%S")
        print "end at", endGuessTime
        startTimeStruct = time.strptime(startGuessTime, "%Y-%m-%d %H:%M:%S")
        endTimeStruct = time.strptime(endGuessTime, "%Y-%m-%d %H:%M:%S")
        timeDiff = time.mktime(endTimeStruct) - time.mktime(startTimeStruct)
        print 'guessing takes ', prettyTimeStr(timeDiff)

        solverLoop(quadSolver) # one more solving step
      
    # end guessing
    sys.stdout.flush()

    # extract key variables
    knownVars = quadSolver.getVars([str(a) for a in pR.gens() if str(a)[0] == "x"])
    bigList = knownVars.split(",") if knownVars else []
    for curItem in bigList:
      var, value = curItem.split("=")
      knownX[var] = value
    print len(bigList), "already known", knownVars


  # check if we are done
  if len(knownX) == 80: break


curStat = quadSolver.stats()
timeSL = int(curStat["time.sl"])
timeEval = int(curStat["time.eval"])
timeDense = int(curStat["time.dense"])
timeSparse = int(curStat["time.sparse"])
timeCPS = int(curStat["time.cps"])

# check if we need to write the system
if args.fullSysName != None:
  print "write full system to file %s"%args.fullSysName
  
  f = open(args.fullSysName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# parameter string was:", sys.argv
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f, "# stats", prettyDictStr(quadSolver.stats())
  print >>f
  print >>f, "# start lin"
  print >>f, quadSolver.getLinSys()
  print >>f
  print >>f, "# start quad"
  print >>f, quadSolver.getQuadSys()
  f.close()

# extract key variables
knownVars = quadSolver.getVars([str(a) for a in pR.gens() if str(a)[0] == "x"])
bigList = knownVars.split(",") if knownVars else []
for curItem in bigList:
  var, value = curItem.split("=")
  knownX[var] = value
print len(bigList), "already known", knownVars

# output the full statistic at the end of the file
print "start at", startTime
endTime = time.strftime("%Y-%m-%d %H:%M:%S")
print "end at", endTime
print "nTriv", nTriv, 
print "ivBlock", args.ivBlock, 
print "Nout=%d, rounds=%d"%(outBits, rounds), 
print "with known x=%d"%len(knownX)

# convert the time
startTimeStruct = time.strptime(startTime, "%Y-%m-%d %H:%M:%S")
endTimeStruct = time.strptime(endTime, "%Y-%m-%d %H:%M:%S")
timeDiff = time.mktime(endTimeStruct) - time.mktime(startTimeStruct)
timeDiff = prettyTimeStr(timeDiff)

print 
print "timeDense timeEval timeSL timeSparse,  / (pretty) timeDense timeEval timeSL timeSparse timeDiff"
print "timeDist", "r%s_o%s_k%s_d%s"%(rounds, outBits, ckey, density), " / ", 
print "%d,%d,%d,%d / %s, %s, %s, %s, %s"%(timeDense, timeEval, timeSL, timeSparse, 
        prettyTimeStr(1.0*timeDense / timeCPS), prettyTimeStr(1.0*timeEval / timeCPS), prettyTimeStr(1.0*timeSL / timeCPS), prettyTimeStr(1.0*timeSparse / timeCPS) , timeDiff)

print 
print "timeDiff repCnt"
print "outCurve", "r%s_o%s"%(rounds,outBits), " / ", 
print "%s,%s"%(timeDiff, repCnt+1)

try:
  # we a done - delete solver
  quadSolver.terminate()
  del quadSolver
  print
  print "============================================================================="
except:
  print "Hello world - I am done (by error)"

