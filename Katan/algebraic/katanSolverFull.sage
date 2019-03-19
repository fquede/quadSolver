#import python useful stuff
import argparse
from sys import *
import os
import time

#import all nesseccary files within this project
from oneKatanFull import *
from eqSolver import *
from utilities import *
from paintMatrix import *

# get external paramters
parser = argparse.ArgumentParser()
# central switches
parser.add_argument('--plainBlock', default="0/256", type=str,
                    help="Start and end for the block of plain text variables. Has the form start/end with start/end being either a number between 0 and 2^80 or a list of iv variables, e.g. i0,1,2. Example: --ivBlock=0/i5. The end is *not* included in the list. Default: 0/i8")
parser.add_argument('--rounds', '-r', default=80, type=int,
                    help="Number of rounds computed from plain. Default: 35. Full Katan uses 254 rounds")
parser.add_argument('--keyGen', '-k', default=1, type=int, 
                    help="1 for standard key, 0 for rnd")
parser.add_argument('--nKat', '-t', default=64, type=int, 
                    help="Number of katan instances generated per solver iteration")

#offline-phase
parser.add_argument('--noOutput', action='store_true', # check
                    help="Suppresses explicit value for output in the system")
parser.add_argument('--writeFile', '-w', default=None, type=str, 
                    help="writes the System to a file (default:None).")

#online-phase
parser.add_argument('--insertOutput', action="store_true", # check
                    help="Add values for the output variables zD_r (D: instance, r: round).")
parser.add_argument('--readFile', '-f', action="append", # check
                    help="Uses a given file as auxiliary input. Multiple files are possible.")
parser.add_argument('--readOrder', '-ro', action="append", # check
                    help="Uses a given file as order.")

#guessing
parser.add_argument('--guess', '-g', default=0, type=int, 
                    help="Guessing up to n variables. 0 for no guessing. Should be less than 40. (default: 0)")
parser.add_argument('--guessRepCnt', '-grc', default=0, type=int,
                    help="0 for guessing at the first solver Iteration. n for guessing at repCnt >= guessrepCnt. There will be just one guessing step! (default: 0)")
parser.add_argument('--guessRead', '-mc', default=None, type=str, 
                    help="sets variables to guess from a file. If None finds variables to guess and writes them to a file (default:None)")
parser.add_argument('--guessWrite', '-gw', default=None, type=str, 
                    help="writes the variables to guess in a file (see below for default)")

parser.add_argument('--toGuess', '-wa', default=[], type=list, 
                    help="determines which variables we are guessing. default: []")

# parameters for quadSolver
parser.add_argument('--density', '-d', default=0.0016, type=float, 
                    help="value between 0 and 1. Gives the minimum density for which we use M4RI")
parser.add_argument('--sl', '-s', default=0, type=int, 
                    help="0 for sparse. 1 for half sparse")
parser.add_argument('--noSL', action="store_true", 
                    help="Switches off sparse linearization (SL) for the solver.")
parser.add_argument('--eval', '-e', default=1, type=int, 
                    help="0 for sparse. 1 for dense")
parser.add_argument('--evalMaxOut', default=None, type=int, 
                    help="Maximal number of monomials when using eval")

#other parameters
parser.add_argument('--paint', '-a', action='store_true',  #check
                    help="activates painting the system matrix (default: False).")


#evalute external parameters
args = parser.parse_args()
rounds = args.rounds
plainBlock = args.plainBlock
ckey = args.keyGen
nKat = args.nKat
guessNum = args.guess
guessFlag = guessNum > 0
density = args.density
slMode = args.sl
evalMode = args.eval
haveFiles = args.readFile != None

#args.toGuess=range(15)+range(19,30)+range(30,40)+range(67,80)
#print args.toGuess,len(args.toGuess)

if haveFiles: 
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file"%curFile

#some internal variables
knownX = dict()
repCntInc = 0
internal = {}
repCnt = 0 

# main solver loop
def solverLoop(quadSolver):
  getOut = False
  evalGoing = True
  keepGoing = 0
  doShortCnt = 1
  while not getOut:
    #main functions
    getOut = True
    keepGoing+=1
    curStats = quadSolver.stats(level=2)
    # we need echelonize now
    if (int(curStats["quadCols"]) == 0) or (int(curStats["quadRows"]) == 0): break # we have a solution
    #curDensity = 1.0
    #curDensity = 0.0
    curDensity = float(curStats["quadNonZeroEntries"])
    curDensity /= float(curStats["quadRows"])*float(curStats["quadCols"])
    if (curDensity < density): quadSolver.echelonizeSparse(); print "S", 
    else: quadSolver.echelonizeDense(); print "D",
    
    # some extra info for the user
    print "(r=%s / c=%s / v=%s)"%(curStats["quadRows"], curStats["quadCols"], curStats["quadVarsInUse"]),
    sys.stdout.flush()
    
    # paint the matrix (if required)
    if args.paint: 
      namePic = "R" + str(rounds) + "_" + "rep" + str(repCnt) + "step" + str(keepGoing) + '_' + '1'
      paintPic(args.paint, quadSolver.getQuadSys(), namePic)

    # eval 
    if evalMode == 0: evalNum = int(quadSolver.fullSparseEval()); print "eS", evalNum,
    else: 
      if args.evalMaxOut == None: evalNum = int(quadSolver.fullEval())
      else: evalNum = int(quadSolver.fullEval(outMax=args.evalMaxOut))  
      print "e", evalNum,
    sys.stdout.flush()
    if evalNum > 0:
      evalGoing = True
      getOut = False
      continue
    else:
      if evalGoing: 
        evalGoing = False
      else: 
        continue

    # sl 
    if not(args.noSL):
      if slMode==0: slNum = quadSolver.sl()
      else: slNum = quadSolver.hsl(curDensity,density)
      if slNum > 2**63: slNum = -(2**64-slNum)
      if slMode==0: print "SL", slNum,
      else: print "hSL", slNum,
      sys.stdout.flush()
      if slNum != 0: getOut = False
  print       
## end of solverLoop

#################################################################################
# tell the world what we are doing
startTime=str(time.strftime("%Y-%m-%d %H:%M:%S"))
print "start at", time.strftime("%Y-%m-%d %H:%M:%S"), "with parameters"
print prettyDictStr(vars(args))
print
sys.stdout.flush()

#################################################################################

# get the order from args.readOrder
filePlainStart = 0; filePlainEnd = 0
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

  filePlainStart = int(min([ var.split('_')[0][1:] for var in varsToQuad ]))
  filePlainEnd = int(max([ var.split('_')[0][1:] for var in varsToQuad ]))
print "start/end of katan instances out of Files", filePlainStart, "/", filePlainEnd

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

#extract the start and end point of wanted katan instances
if args.plainBlock != None: 
  plainStart,plainEnd = args.plainBlock.split("/")
  plainStart = int(plainStart); plainEnd = int(plainEnd)

# create a set of all Katan instances that should be generated
if filePlainEnd!=0: filePlainEnd+=1
knownKatan = set(range(filePlainStart,filePlainEnd))
wantedKatan = set(range(plainStart, plainEnd))

repPlainSet = sorted(wantedKatan.difference(knownKatan))
# check how many instances we have in our current system (including files)
curEnd = max(filePlainEnd,plainEnd)
plainVarNum = 0
while 2^plainVarNum <= curEnd: plainVarNum += 1
#decide when to guess
if args.guessRepCnt == 0: args.guessRepCnt = curEnd

# key setup
if ckey==1: KE = [0 if ((a % 3) == 0) or ((a % 7) == 0) else 1 for a in xrange(80)] #standard key
#if ckey==1: KE = [GF(2)(1) for _ in xrange(80)] #1 key
else: KE = [GF(2).random_element() for _ in xrange(80)]# random key

#setting up the BooleanPolynomialRing 
keyVars = [ "x%02d"%i for i in xrange(80) ]
#roundKeyVars = [ 'k%03d'%i for i in xrange(2*rounds)]
stateVars = flatten([ ["aD_%04d"%i,"bD_%04d"%i] for i in xrange(rounds) ])
helpVars = [ "gD_%02d"%i for i in xrange(12) ]
helpVars += [ "hD_%02d"%i for i in xrange(18) ]
varNames = ['cD_%02d'%j for j in range(32)] + stateVars + helpVars + [ "pD_%02d"%i for i in xrange(32) ] + keyVars
varNames.reverse()
pR = BooleanPolynomialRing(len(varNames), varNames, order='degneglex')

#reverse the keyVars for later use
#keyVars.reverse()
#roundKeyVars.reverse()

print "start=%d/end=%d"%(plainStart, plainEnd),
print "number of instances to generate: ", len(repPlainSet)

# plain & key variables in a handy dictionary
valueDict = {}
for i in range(80):
  valueDict[pR("x%02d"%i)] = KE[i]
print "the key is: ", prettyDictStr(valueDict)
for i in range(plainVarNum):
  valueDict[pR("pD_%02d"%i)] = 0 # for the time being - will be changed later

# create "golden instance"
katanEqSysPlain = oneKatan32enc(rounds, pR, plainVarNum,valueDict,args.toGuess)
katanEqSysCipher = oneKatan32dec(rounds, pR, plainVarNum,valueDict,args.toGuess)

# move the system of equations to an array of strings
katanEqSysPlain = [ str(a).replace(" ", "") for a in katanEqSysPlain ]
#katanEqSysPlain = [  ]
katanEqSysCipher = [ str(a).replace(" ", "") for a in katanEqSysCipher ]
#for f in katanEqSysCipher:
#  print f
#katanEqSysCipher = [  ]

# statistic
fullStat = {}
fullStat["x"] = -1
fullStat["vars"] = -1
fullStat["eqs"] = -1
fullStat["repCnt"] = -1
fullStat["delta"] = -1
fullStat["mons"] = -1
fullStat["ratio"] = -1
monAbs = 0; monMax = -1

# set up solver
quadSolver = eqSolver()

#round-key and key variables have to be put first into the solver
#quadSolver.addVariables(roundKeyVars)
#quadSolver.addVariables(keyVars)

# list of new variables for every instance
varsToQuad.update(keyVars)
varsToQuad.update([ "c%04d_%02d"%(i,j) for j in range(32) for i in range(plainStart, plainEnd) ])
varsToQuad.update([ "p%04d_%02d"%(i,j) for j in range(32) for i in range(plainStart, plainEnd) ])
for v in repPlainSet:
  varsToQuad.update([ y.replace("D",'%04d'%v) for y in stateVars ])
varsToQuad = sorted(varsToQuad, key=lambda x: varname2key(rounds, x))
quadSolver.addVariables(varsToQuad) # define variables for the solver

#adding key to quadSolver; debug purposes only!
#quadSolver.addEquations(['x%02d+%d'%(i,KE[i]) for i in range(80)])

#storing the order in a file
if args.writeFile:
  orderFileName = args.writeFile+'.order'
  f = open(orderFileName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f
  print >>f, "# start order"
  for curVar in varsToQuad:
    print >>f, curVar
  f.close()

# add output values
if args.insertOutput:
  byInstance = {}
  varDict = { "x%02d"%i: KE[i] for i in range(80) }
  for instanceNum in knownKatan:
    byInstance[instanceNum] = [plainVarNum, varDict, args.rounds]
  print "add output"
  sys.stdout.flush()
  # generate the corresponding output bit in parallel
  for X,Y in katanOutVals(byInstance.items()):
    quadSolver.addEquations(Y)
  print "--- added values for output"

#add equations from Files to the quadSolver
if args.readFile:
  quadSolver.addEquations(linEqs)
  quadSolver.addEquations(quadEqs)

#adding equations for the known vars in the quadSolver
#quadSolver.addEquations(['x%02d+'%i + str(valueDict[pR('x%02d'%i)]) for i in args.toGuess])

#######################
# delete variables
del varsToQuad, linEqs, quadEqs

# set the maximal number of rounds to terminate the algorithm
katCount = 5

if repPlainSet == []: repPlainSet.append('oI') #make sure there is an iteration
for repCnt in repPlainSet:
  fullStat["repCnt"] = repCnt
  print "repCnt", repCnt, 

  # specialize the instance
  if repCnt != 'oI':
    curInstance = specialInstance(repCnt, katanEqSysPlain, katanEqSysCipher, rounds, valueDict, plainVarNum, args.noOutput)
    quadSolver.addEquations(curInstance)

  if repCnt == 'oI': repCnt = args.guessRepCnt # make sure we make a full iteration
  
  # start solving - if we have enough instances (or are in the last round anyway)
  if ((repCnt % katCount) == (katCount-1)) or (repCnt >= max(repPlainSet)) or (repCnt >= len(repPlainSet)):
    print
    print "rep", repCnt, "time", time.strftime("%Y-%m-%d %H:%M:%S")
    print "   stats before are", prettyDictStr(quadSolver.stats(level=2))
    sys.stdout.flush()

    #setting katCount (big step / baby step)
    if (katCount == nKat) and (nKat>3): katCount += floor(nKat/4)
    else: katCount = nKat

    # start solving
    if not(guessFlag) or (repCnt <= args.guessRepCnt - 1): solverLoop(quadSolver)

    # check if we need to write the system
    if args.writeFile != None:
      print "write full system to file %s"%args.writeFile
      
      f = open(args.writeFile+str(repCnt), "w")
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
    if guessFlag and (repCnt >= args.guessRepCnt-1):
      guessFlag = False  # only guess once
      startGuessTime=str(time.strftime("%Y-%m-%d %H:%M:%S"))
      print 'start guessing at ',startGuessTime

      #getting the quadratic systems out of the Solver
      quadSys = quadSolver.getQuadSys().replace(" ", "").split("\n")
      quadSys.pop() # get rid of the zero
      #getting the variables to guess
      if args.guessRead:
        print "reading vars to guess out of File: ", args.guessRead
        f = open(args.guessRead,'r')
        control_string = f.readline()
        print control_string
        control_string = control_string.split(',')
        assert (control_string[0][2:].strip() == str(rounds) and control_string[1][8:].strip() <= str(args.guessRepCnt)), \
                "variables to guess are not for the current problem"
        varList = eval(f.read())
        f.close()
      else:
        # double-counting the variables
        maxVars = []
        #counting per poly
        for poly in quadSys:
          allVars = flatten([term.split('*') for term in poly.split('+')])

          cntDict = { a:0 for a in allVars }
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
        f.write('R='+str(rounds)+',gRepCnt='+str(args.guessRepCnt)+'\n')
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
        #guessedVars = []
        #for curVar in varToGuess:
        #  guessedVars.append(curVar + '+%d'%(GF(2).random_element()))
        byInstance = {}
        guessedVars = []
        varDict = { "x%02d"%i:KE[i] for i in range(80) }
        for curVar in varToGuess:
          if curVar.startswith('x'):
            guessedVars.append(curVar + '+%d'%varDict[curVar])
            continue
          if curVar.startswith('k'):
            key_lfsr = [ GF(2)(KE[i]) for i in xrange(80) ]
            for _ in xrange(int(curVar[1:])+1):
              round_key = key_lfsr[0]
              new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]
              key_lfsr[:79] = key_lfsr[1:]
              key_lfsr[79] = new_bit
              
            guessedVars.append(curVar + '+%d'%round_key)
            continue
          instance = curVar.split("_")[0][1:]
          if not(byInstance.has_key(instance)): byInstance[instance] = [plainVarNum,varDict,set()]
          byInstance[instance][2].add(curVar)

        # generate the corresponding data in parallel
        for X,Y in katanData(byInstance.items()):
          guessedVars += Y

        #adding the guessed values to the Solver
        print "there are %d key variables from"%len(filter(lambda v: v[0] == "x", guessedVars)), guessedVars 
        quadSolver.addEquations([s.replace("=", "+") for s in guessedVars])     

        #time used for guessing 
        endGuessTime = time.strftime("%Y-%m-%d %H:%M:%S")
        print "we are done guessing at", endGuessTime
        startTimeStruct = time.strptime(startGuessTime, "%Y-%m-%d %H:%M:%S")
        endTimeStruct = time.strptime(endGuessTime, "%Y-%m-%d %H:%M:%S")
        timeDiff = time.mktime(endTimeStruct) - time.mktime(startTimeStruct)
        print 'guessing takes ', prettyTimeStr(timeDiff)

        solverLoop(quadSolver) # one more solving step
      
    # end guessing
    sys.stdout.flush()

    # extract key variables
    knownVars = quadSolver.getVars(['x%02d'%i for i in xrange(80) ])
    bigList = knownVars.split(",") if knownVars else []
    for curItem in bigList:
      var, value = curItem.split("=")
      knownX[var] = value
    print len(bigList), "already known", knownVars
  sys.stdout.flush()

  # check if we are done
  if len(knownX) >= 80: break


curStat = quadSolver.stats()
timeSL = int(curStat["time.sl"])
timeEval = int(curStat["time.eval"])
timeDense = int(curStat["time.dense"])
timeSparse = int(curStat["time.sparse"])
timeCPS = int(curStat["time.cps"])

# extract key variables
knownVars = quadSolver.getVars(['x%02d'%i for i in xrange(80)])
bigList = knownVars.split(",") if knownVars else []
for curItem in bigList:
  var, value = curItem.split("=")
  knownX[var] = value
print len(bigList), "already known", knownVars


# we a done - delete solver
quadSolver.terminate()

# output the full statistic at the end of the file
print "start at", startTime
endTime = time.strftime("%Y-%m-%d %H:%M:%S")
print "end at", endTime
print "nKat", nKat
print "plainBlock", plainBlock
print "stat collected:", fullStat
print "rounds=%d"%rounds, 
print "with known x=%d from"%len(knownX), 
print "#vars=%d, #eqs=%d, #delta=%d, #mons=%s, ratio=%.3f"%(fullStat["vars"], fullStat["eqs"], fullStat["delta"], fullStat["mons"], fullStat["ratio"])

# convert the time
startTimeStruct = time.strptime(startTime, "%Y-%m-%d %H:%M:%S")
endTimeStruct = time.strptime(endTime, "%Y-%m-%d %H:%M:%S")
timeDiff = time.mktime(endTimeStruct) - time.mktime(startTimeStruct)
timeDiff = prettyTimeStr(timeDiff)

print
print "repCnt, timeDiff"
print "plainJob", "r%s_k%s"%(rounds, ckey), " / ", 
print "%s,%s"%(repCnt+1, timeDiff)

print 
print "timeDense timeEval timeSL timeSparse,  / (pretty) timeDense timeEval timeSL timeSparse timeDiff"
print "timeDist", "r%s_k%s_d%s"%(rounds, ckey, density), " / ", 
print "%d,%d,%d,%d / %s, %s, %s, %s, %s"%(timeDense, timeEval, timeSL, timeSparse, 
        prettyTimeStr(1.0*timeDense / timeCPS), prettyTimeStr(1.0*timeEval / timeCPS), prettyTimeStr(1.0*timeSL / timeCPS), prettyTimeStr(1.0*timeSparse / timeCPS) , timeDiff)

try:
  print
  print "============================================================================="
except:
  print "Hello world - I am so close to being a good solver"

