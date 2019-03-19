#import python useful stuff
import argparse
from sys import *
import os
import time

#import all nesseccary files within this project
from twoKatan import *
from eqSolver import *
from utilities import *
from paintMatrix import *

# get external paramters
parser = argparse.ArgumentParser()
# central switches
parser.add_argument('--plainVars', '-pv', default=8, type=int,
                    help="Number of variables used in the plaintext")
parser.add_argument('--rounds', '-r', default=254, type=int,
                    help="Number of rounds computed for Katan. Full Katan uses 254 rounds")
parser.add_argument('--faultRound', '-fR', default=246, type=int,
                    help='Number of the round in which the fault is inserted')
parser.add_argument('--rounddelta', '-Rd', default=8, type=int,
                    help='Round delta')
parser.add_argument('--keyGen', '-k', default=0, type=int, 
                    help="1 for standard key, 0 for rnd")
parser.add_argument('--maxFault', '-m', default=20, type=int, 
                    help="Number of maximal Fault injections")
parser.add_argument('--faultModel', '-fm', default='3f', type=str, 
                    help="faultModel which is used. 1f for a one fault model which changing position and 3f for random vector which is set after the first occurence.")
parser.add_argument('--faultVec', '-fv', action="append", # check
                    help="Uses a given file as order.")

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
parser.add_argument('--paint', '-a', action='store_true',  
                    help="activates painting the system matrix (default: False).")
parser.add_argument('--writeSys', '-w', default=None, type=str,
                    help="the name of the File where you can find the resulting system (default: None)")


#evalute external parameters
args = parser.parse_args()
rounds = args.rounds
numPlainVars = args.plainVars
ckey = args.keyGen
density = args.density
slMode = args.sl
evalMode = args.eval
faultRound=args.faultRound

#checking the fault model
assert args.faultModel=='1f' or args.faultModel=='3f'

#some internal variables
knownX = dict()
repCntInc = 0
internal = {}
repCnt = 0 

# main solver loop
def solverLoop(loopSolver):
  getOut = False
  evalGoing = True
  keepGoing = 0
  doShortCnt = 1
  while not getOut:
    #main functions
    getOut = True
    keepGoing+=1
    curStats = loopSolver.stats(level=2)
    # we need echelonize now
    if (int(curStats["quadCols"]) == 0) or (int(curStats["quadRows"]) == 0): break # we have a solution
    #curDensity = 1.0
    #curDensity = 0.0
    curDensity = float(curStats["quadNonZeroEntries"])
    curDensity /= float(curStats["quadRows"])*float(curStats["quadCols"])
    if (curDensity < density): loopSolver.echelonizeSparse(); print "S", 
    else: loopSolver.echelonizeDense(); print "D",
    
    # some extra info for the user
    print "(r=%s / c=%s / v=%s)"%(curStats["quadRows"], curStats["quadCols"], curStats["quadVarsInUse"]),
    sys.stdout.flush()
    
    # paint the matrix (if required)
    if args.paint: 
      namePic = "R" + str(rounds) + "_" + "fR" + str(args.faultRound) + "step" + str(keepGoing)
      paintPic(args.paint, loopSolver.getQuadSys(), namePic)

    # eval 
    if evalMode == 0: evalNum = int(loopSolver.fullSparseEval()); print "eS", evalNum,
    else: 
      if args.evalMaxOut == None: evalNum = int(loopSolver.fullEval())
      else: evalNum = int(loopSolver.fullEval(outMax=args.evalMaxOut))  
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
      if slMode==0: slNum = loopSolver.sl()
      else: slNum = loopSolver.hsl(curDensity,density)
      if slNum > 2**63: slNum = -(2**64-slNum)
      if slMode==0: print "SL", slNum,
      else: print "hSL", slNum,
      sys.stdout.flush()
      if slNum != 0: getOut = False
  print       
## end of solverLoop
#################################################################################

def multSolver(origSolver, faultPos, faultRound, numFault, newFault, faultVecs):
  #changing the fault round
  faultRound-=args.rounddelta
  if faultRound<=182: 
    #oldFaultPos=faultPos
    faultRound = 245
    #faultPos=ZZ.random_element(32)
    print 'new position: ', faultPos

  numFault+=1
  fullStat["numFault"] = numFault
  
  #abort numFault
  if numFault > args.maxFault: return False, origSolver, numFault, False, faultVecs, faultPos, faultRound
  print numFault, 'faults done. This fault in round', faultRound
  
  
  if newFault:
    faultVec=[] 
    Fbit1=0;Fbit2=0;Fbit3=0
    while (Fbit1==0 and Fbit2==0 and Fbit3==0) or (faultPos==0 and Fbit2==0 and Fbit3==0) or (faultPos==31 and Fbit1==0 and Fbit2==0):
      Fbit1,Fbit2,Fbit3 = [GF(2).random_element() for _ in range(3)]
    if fault_pos==0:
      if Fbit2==1: faultVec.append(fault_pos)
      if Fbit3==1: faultVec.append(fault_pos + 1)
    elif fault_pos==31:
      if Fbit1==1: faultVec.append(fault_pos - 1)
      if Fbit2==1: faultVec.append(fault_pos)
    else:
      if Fbit1==1: faultVec.append(fault_pos - 1)
      if Fbit2==1: faultVec.append(fault_pos)
      if Fbit3==1: faultVec.append(fault_pos + 1)
    print 'new fault vector is:', faultVec
    faultVecs.append(faultVec)
  else:
    faultVec = faultVecs[numFault-2]

  #creating all vectors
  if faultPos==0:
    tryVecs=[[0],[1],[0,1]]
  elif faultPos==31:
    tryVecs=[[30],[31],[30,31]]
  else: tryVecs=[[faultPos-1],[faultPos],[faultPos+1],[faultPos-1, faultPos],[faultPos, faultPos+1],[faultPos-1, faultPos+1],[faultPos-1, faultPos, faultPos+1]]
  
  
  for tryVec in tryVecs:
    print 'trying vector: ', tryVec
    #copying the original solver
    trySolver = eqSolver()
    trySolver.addVariables(varsToQuad) # define variables for the solver
    trySolver.addEquations(origSolver.getLinSys().replace(" ","").split('\n'))
    trySolver.addEquations(origSolver.getQuadSys().replace(" ","").split('\n'))

    #generating the system
    curSystem = twoKatanVec(pR, KE, plain, rounds, faultRound, faultVec, tryVec)
    # move the system of equations to an array of strings
    curSystem = [ str(a).replace(" ", "") for a in curSystem ]
    curSystem = [ s.replace("F", '%04d'%numFault) for s in curSystem]
    

    #adding the system to the new Solver
    trySolver.addEquations(curSystem)
    #start solving
    solverLoop(trySolver)
    
    #write the system
    if args.writeSys != None:
      print "write full system to file %s"%(args.writeSys+str(faultRound)+'_'+str(numFault))
      
      f = open(args.writeSys+str(faultRound)+'_'+str(numFault), "w")
      print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
      print >>f, "# parameter string was:", sys.argv
      print >>f, "# all parameters were:", prettyDictStr(vars(args))
      print >>f, "# stats", prettyDictStr(trySolver.stats())
      print >>f
      print >>f, "# start lin"
      print >>f, trySolver.getLinSys()
      print >>f
      print >>f, "# start quad"
      print >>f, trySolver.getQuadSys()
      f.close()

    # some statistics
    curStats = trySolver.stats(level=2)
    print "stats after are", prettyDictStr(curStats)
    sys.stdout.flush()

    # extract key variables
    knownX={}
    knownVars = trySolver.getVars(['x%02d'%i for i in xrange(80) ])
    bigList = knownVars.split(",") if knownVars else []
    for curItem in bigList:
      var, value = curItem.split("=")
      knownX[var] = value
    print len(bigList), "already known", knownVars
    
    # check if we are done
    if len(knownVars) != 0:
      if knownVars[0]=='!': continue #incorrect solution
    if len(knownX) == 80: return True, trySolver, numFault, False, faultVecs, faultPos, faultRound # correct solution

    #not done yet
    solFlag, trySolver, numFault, newFault, faultVecs, faultPos, faultRound = multSolver(trySolver,faultPos,faultRound,numFault,True,faultVecs)
    if solFlag == True: return True, trySolver, numFault, False, faultVecs, faultPos, faultRound
  #all wrong
  if faultRound==246: faultRound=190 #; faultPos=oldFaultPos
  else: faultRound+=args.rounddelta
  return False, origSolver, numFault-1, False, faultVecs, faultPos, faultRound
    
    

#################################################################################
# tell the world what we are doing
startTime=str(time.strftime("%Y-%m-%d %H:%M:%S"))
print "start at", time.strftime("%Y-%m-%d %H:%M:%S"), "with parameters"
print prettyDictStr(vars(args))
print
sys.stdout.flush()

#################################################################################

# key and plain text setup
if ckey==1: KE = [0 if ((a % 3) == 0) or ((a % 7) == 0) else 1 for a in xrange(80)] #standard key
#if ckey==1: KE = [GF(2)(1) for _ in xrange(80)] #1 key
else: KE = [GF(2).random_element() for _ in xrange(80)]# random key
plain = [GF(2).random_element() for _ in range(numPlainVars)] + [GF(2)(0) for _ in range(32-numPlainVars)]

#setting up the BooleanPolynomialRing 
keyVars = [ "x%02d"%i for i in xrange(80) ]
roundKeyVars = [ 'k%03d'%i for i in xrange(2*rounds)]
stateVars = flatten([ ["aD_%04d"%i,"bD_%04d"%i] for i in xrange(rounds) ])
stateVarsFaulty = flatten([ ["aF_%04d"%i,"bF_%04d"%i] for i in xrange(rounds) ])
varNames = ['zF_%02d'%j for j in range(32)] + ['zD_%02d'%j for j in range(32)] + stateVarsFaulty + stateVars + [ "pD_%02d"%i for i in xrange(32) ] + keyVars + roundKeyVars + ['y%02d'%i for i in range(32)]
varNames.reverse()
pR = BooleanPolynomialRing(len(varNames), varNames, order='degneglex')
  
# plain & key variables in a handy dictionary
valueDict = {}
for i in range(80):
  valueDict["x%02d"%i] = KE[i]
print "the key is: ", prettyDictStr(valueDict)
for i in range(numPlainVars):
  valueDict["pD_%02d"%i] = plain[i] 

# statistic
fullStat = {}
fullStat["x"] = -1
fullStat["vars"] = -1
fullStat["eqs"] = -1
fullStat["numFault"] = -1
fullStat["delta"] = -1
fullStat["mons"] = -1
fullStat["ratio"] = -1
monAbs = 0; monMax = -1

# set up solver
quadSolver = eqSolver()

# list of new variables for every instance
global varsToQuad
varsToQuad=set(roundKeyVars)
varsToQuad.update(keyVars)
varsToQuad.update([ "zD_%02d"%j for j in range(32) ] + ['y%02d'%i for i in range(32)])
varsToQuad.update(stateVars)
for v in range(args.maxFault):
  varsToQuad.update( [ y.replace("F",'%04d'%v) for y in [ "zF_%02d"%j for j in range(32) ]])
  varsToQuad.update([ y.replace("F",'%04d'%v) for y in stateVarsFaulty ])
varsToQuad = sorted(varsToQuad, key=varname2key)
quadSolver.addVariables(varsToQuad) # define variables for the solver

#debugging only
#quadSolver.addEquations( ['x%02d+'%i+str(KE[i]) for i in range(80) ] )

list_to_guess=[4,7,2,5,6,9,3,0,14,15,12,10,19,13,18,11,20,23,1,8,31,22,21,38,28,26,24,17,30,16,39,27,47,43,29,34,35,37]
#11,0,1,9,3,5,7,13,2,6,15,14,4,16,8,17,190, 24, 12, 29, 18, 21, 23, 20, 30, 19, 33, 32, 25, 28, 27, 34, 22, 26, 35, 43]
print len(list_to_guess)
quadSolver.addEquations( ['x%02d+'%i+str(KE[i]) for i in list_to_guess ] )


#start injecting Faults and solving the System
if args.faultModel=='3f': faultRound=args.faultRound

numFault=0
faultVec = args.faultVec
fault_pos = ZZ.random_element(32)
print 'position is: ', fault_pos
if faultVec==[]: 
  faultVec=None
if faultVec:
  for i in range(len(faultVec)):
    faultVec[i] = int(faultVec[i])
print faultVec

while numFault < args.maxFault:
  numFault+=1
  fullStat["numFault"] = numFault
  print numFault, 'faults done. This fault in round', faultRound

  #generating the system with one Fault ########## dynamic faultRound?!?
  if args.faultModel=='1f':
    curSystem=twoKatan(pR, KE, plain, rounds, faultRound)
  else: 
    if faultRound<241 and args.faultModel=='3f':
      solFlag, quadSolver, numFault, newFault, faultVecs, faultPos, faultRound = multSolver(quadSolver, fault_pos, faultRound+args.rounddelta, numFault-1, True, [])
      break
    else: 
      faultVec=[] 
      Fbit1=0;Fbit2=0;Fbit3=0
      while (Fbit1==0 and Fbit2==0 and Fbit3==0) or (fault_pos==0 and Fbit2==0 and Fbit3==0) or (fault_pos==31 and Fbit1==0 and Fbit2==0):
        Fbit1,Fbit2,Fbit3 = [GF(2).random_element() for _ in range(3)]
      if fault_pos==0:
        if Fbit2==1: faultVec.append(fault_pos)
        if Fbit3==1: faultVec.append(fault_pos + 1)
      elif fault_pos==31:
        if Fbit1==1: faultVec.append(fault_pos - 1)
        if Fbit2==1: faultVec.append(fault_pos)
      else:
        if Fbit1==1: faultVec.append(fault_pos - 1)
        if Fbit2==1: faultVec.append(fault_pos)
        if Fbit3==1: faultVec.append(fault_pos + 1)
      print 'the fault vector is:', faultVec
      curSystem = twoKatanVec(pR, KE, plain, rounds, faultRound, faultVec, None)
    
  #changing fault round
  if args.faultModel=='1f':
    faultRound-=8
    if numFault%9==0: faultRound = args.faultRound
  else:
    faultRound-=18
    if faultRound<=182: 
      faultRound = 246
      fault_pos=ZZ.random_element(32)
      print 'new position: ', fault_pos
  
  # move the system of equations to an array of strings
  curSystem = [ str(a).replace(" ", "") for a in curSystem ]
  curSystem = [ s.replace("F", '%04d'%numFault) for s in curSystem]
  #print curSystem
  #adding the system to the Solver
  quadSolver.addEquations(curSystem)

  #start solving
  solverLoop(quadSolver)

  #write the system
  if args.writeSys != None:
    print "write full system to file %s"%(args.writeSys+str(numFault))
    
    f = open(args.writeSys+str(numFault), "w")
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

  #extract roundkeys
  knownVarsK = quadSolver.getVars(['k%03d'%i for i in xrange(508) ])
  bigListK = knownVarsK.split(",") if knownVarsK else []
  #for curItem in bigListK:
  #  var, value = curItem.split("=")
  #  knownK[var] = value
  print len(bigListK), "already known", knownVarsK
  sys.stdout.flush()
  if len(bigListK) >= 80: solverLoop(quadSolver)

  # extract key variables
  knownVars = quadSolver.getVars(['x%02d'%i for i in xrange(80) ])
  bigList = knownVars.split(",") if knownVars else []
  for curItem in bigList:
    var, value = curItem.split("=")
    knownX[var] = value
  print len(bigList), "already known", knownVars
  if len(knownVars) != 0:
    if knownVars[0]=='!': break

  # check if we are done
  if len(knownX) >= 80: break

# extract key variables
knownX={}
knownVars = quadSolver.getVars(['x%02d'%i for i in xrange(80) ])
bigList = knownVars.split(",") if knownVars else []
for curItem in bigList:
  var, value = curItem.split("=")
  knownX[var] = value
print len(bigList), "already known", knownVars

# extract char variables
#knownX={}
#knownVars = quadSolver.getVars(['y%02d'%i for i in xrange(32) ])
#print knownVars
#bigList = knownVars.split(",") if knownVars else []
#print bigList
#
#charString = ''
#notEval = 0
#for i in range(32):
#  if 'y%02d'%i == bigList[i-notEval][0:3]: charString+=bigList[i-notEval][-1]+' '
#  else: charString+= '- '; notEval+=1
#print charString

#for curItem in bigList:
#  var, value = curItem.split("=")
#  knownX[var] = value
#print len(bigList), "already known", knownVars

curStat = quadSolver.stats()
timeSL = int(curStat["time.sl"])
timeEval = int(curStat["time.eval"])
timeDense = int(curStat["time.dense"])
timeSparse = int(curStat["time.sparse"])
timeCPS = int(curStat["time.cps"])

# extract key variables
#knownVars = quadSolver.getVars(['x%02d'%i for i in xrange(80)])
#bigList = knownVars.split(",") if knownVars else []
#for curItem in bigList:
#  var, value = curItem.split("=")
#  knownX[var] = value
#print len(bigList), "already known", knownVars


# we a done - delete solver
quadSolver.terminate()

# output the full statistic at the end of the file
print "start at", startTime
endTime = time.strftime("%Y-%m-%d %H:%M:%S")
print "end at", endTime
print "maximal number of Faults", args.maxFault
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
print "numFault, timeDiff"
print "faulty", "r%s_k%s"%(rounds, ckey), " / ", 
print "%s,%s"%(numFault, timeDiff)

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

