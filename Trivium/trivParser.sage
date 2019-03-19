# frontEnd of the solver for Trivium
# Use via command line (options see below)

#import python useful stuff
import argparse
from sys import *
import os
import time

# import / attach all necessary files within the project
import configTrivSolver

from trivAutomat import *

########################################################################
# get external paramters
parser = argparse.ArgumentParser()

#instances
parser.add_argument('--instances', action='store_true',
                    help="Returning a system for the trivium instances according to the rounds, outBits and ivBlock parameter")
parser.add_argument('--rounds', '-r', default=560, type=int,
                    help="Number of rounds computed before producing output. Default: 560. Full Trivium uses 1152 rounds")
parser.add_argument('--outBits', '-o', default=1, type=int,
                    help="Number of output rounds used in the equations. Default: 1")
parser.add_argument('--writeName', default=None, type=str, 
                    help="If set, writes the full system to this file. Overwrites another file of the same name that may already exist.")
parser.add_argument('--ivBlock', default="0/256", type=str,
                    help="Start and end for the block of IV variables. Has the form start/end with start/end being either a number between 0 and 2^80. As for range, the end is *not* included in the list. Default: 0/256")

#merge #writeName
parser.add_argument('--merge', action='store_true',
                    help="Merging two generated systems with the same round into one.")
parser.add_argument('--readFile', '-f', action="append", 
                    help="Uses a given file as auxiliary input. Multiple files are possible.")
parser.add_argument('--cubeLen', default=150, type=int,
                    help="Maximal number of terms in a cube. Default: 150")

#data collection #writeName, rounds, outBits, ivBlock
parser.add_argument('--data', action='store_true',
                    help="Producing the data (output and guessing variables) from a given or random key.")
parser.add_argument('--guessTest', action='store_true',
                    help="testing the fail time rate of the solver")
parser.add_argument('--minCover', '-mc', default=None, type=str, 
                    help="sets variables to guess from a file. If None finds variables to guess and writes them to a file (default:None)")
parser.add_argument('--keyGen', '-k', default=1, type=int, 
                    help="1 for standard key, 0 for rnd")
parser.add_argument('--guess', '-g', default=0, type=int, 
                    help="Guessing up to n variables. 0 for no guessing. Should be less than 40. (default: 0)")

#export #writeName, readFile
parser.add_argument('--export', action='store_true',
                    help="Producing the data (output and guessing variables) from a given or random key.")
parser.add_argument('--guessWrite', '-gw', default=None, type=str,
                    help="writes the variables to guess in a file (see below for default)")
parser.add_argument('--readOrder', '-ro', default=None, type=str,
                    help="writes the variables to guess in a file (see below for default)")

#solve, #readFile, readOrder, guessTest
parser.add_argument('--solve', action='store_true',
                    help="Solving a given system")
parser.add_argument('--readData', '-rd', default=None, type=str,
                    help="reads the data out of the data collection section")

#check consistency, #readFile, readOrder, readData
parser.add_argument('--checkConsistency', action='store_true',
                    help="Checks the consistency of the system")

# strengthen, #readFile, readOrder, writeName
parser.add_argument('--strengthen', action='store_true',
                    help="Given a system, applies Sparse Target Linearization (STL) and a Genetic Algorithm (GA).")
parser.add_argument('--insertCubes', action="append", 
                    help="Adds additional cubes equations. Multiple files are possible")
# Mutants / stl (sparse target linearization)
parser.add_argument('--stlMaxTime', default=3*3600, type=int, 
                    help="Maximal running time of the stl subsystem in seconds. Default: 3h.")
parser.add_argument('--stlMaxMons', default=10^5, type=int, 
                    help="Maximal number of monomials considered in one run of the outer loop. Default: 10^5.")
parser.add_argument('--stlMaxNewMons', default=250, type=int, 
                    help="Maximal number of new target monomials. Default: 250.")
parser.add_argument('--stlMaxSingleLen', default=30, type=int, 
                    help="Maximal length of polynomials with exactly one degree 2 term. Default: 30.")
parser.add_argument('--stlMaxMultiLen', default=10, type=int, 
                    help="Maximal length of polynomials with at least two degree 2 terms. Default: 10.")
parser.add_argument('--stlHard', action='store_true',
                    help="Allows new degree 2 monomials in the stl step.")
# gps: genetic polynomial search
parser.add_argument('--gpsMaxCopy', default=10^5, type=int, 
                    help="Maximal number of polynomials that are used as input for the GA. Default: 10^5.")
parser.add_argument('--gpsRepNum', default=25, type=int, 
                    help="Maximal number of (internal) repetitions of the GA. Default: 25.")
parser.add_argument('--gpsMaxPolyLen', default=30, type=int, 
                    help="Maximal number of terms for the polynomials in the GA. Default: 30.")
parser.add_argument('--gpsResSize', default=10^5, type=int, 
                    help="Maximal polynomials considered in the GA. Default: 10^5.")
parser.add_argument('--gpsInsertMem', default=4, type=float, 
                    help="Maximal amount of extra memory for inserting polynomials in GB. Default: 4.")
# cnf: conjugated normal form (preprocessor)
parser.add_argument('--cnfOn', action='store_true',
                    help="Activates the anf <-> cnf converstion routines.")
parser.add_argument('--cnfMaxRead', default=5, type=int, 
                    help="Maximal number of variables per clause that are transferred back to ANF. Default: 5.")
parser.add_argument('--cnfMaxWrite', default=20, type=int, 
                    help="Maximal number of terms per polynomial that are transferred to CNF. Default: 20.")

# parameters for quadSolver
parser.add_argument('--evalMaxOut', default=None, type=int, 
                    help="Maximal number of monomials when using eval")
parser.add_argument('--density', '-d', default=0.0016, type=float, 
                    help="value between 0 and 1. Gives the minimum density for which we use M4RI")
parser.add_argument('--sl', '-s', default=0, type=int, 
                    help="0 for sparse. 1 for half sparse. -1 disables SL")
parser.add_argument('--evalSparse', '-es', action="store_true", 
                    help="enables sparse evaluation.")

########################################################################


#evalute external parameters
args = parser.parse_args()
guessFlag = args.guess > 0
haveFiles = args.readFile != None

# check if all necessary files are present
if haveFiles: 
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file (readFile)"%curFile
if args.insertCubes != None: 
  for curFile in args.insertCubes:
    assert os.path.isfile(curFile), "%s is not a file (insetCubes)"%curFile

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

# make sure there is a file name
if not(args.solve):
  if args.writeName == None:
    args.writeName = 'r' + str(args.rounds) + '_o' + str(args.outBits) + '_b' + str(args.ivBlock).replace("/", "_") + '.txt'
    print 'No file name is given. Using ' + args.writeName
  else:
    print "write full system to file %s"%args.writeName
# issue a warning if necessary
if args.writeName != None:
  if os.path.isfile(args.writeName): print "!!! %s does already exist !!! Will be overwritten !!!"%args.writeName

# start with triv solving
if args.instances: instances(args)

if args.merge:
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file"%curFile
  merge(args)

if args.strengthen:
  # apply mutants and genetic algorithms to derive a stronger system
  strengthen(args)

if args.export:
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file"%curFile
  export(args)

if args.data:
  for curFile in args.readFile:
    assert os.path.isfile(curFile), "%s is not a file"%curFile
  data_collection(args)

# just solving the system
if args.solve: solveSys(args)

if args.checkConsistency:
  # just solving the system
  check_consistency(args)

# output the full statistic at the end of the file
print "start at", startTime
endTime = time.strftime("%Y-%m-%d %H:%M:%S")
# convert the time
startTimeStruct = time.strptime(startTime, "%Y-%m-%d %H:%M:%S")
endTimeStruct = time.strptime(endTime, "%Y-%m-%d %H:%M:%S")
timeDiff = time.mktime(endTimeStruct) - time.mktime(startTimeStruct)
timeDiff = prettyTimeStr(timeDiff)
print "end at", endTime
print "time difference ", timeDiff
print "ivBlock", args.ivBlock, 
print "Nout=%d, rounds=%d"%(args.outBits, args.rounds)
exit(0)


# verbosity set
if args.verbose != None: 
  verbStr = ",".join(args.verbose)
  verbositySet = set(verbStr.split(","))
else: verbositySet = set()


