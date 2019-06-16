import multiprocessing as mp
import os
from time import time

solutionDir = 'C:\\DOPT-2019-Practice'
hwName = 'hw2'
taskName = 'tsp'
buildPath = os.path.join(solutionDir, hwName,'x64\\Release\\')

solverPath = os.path.join(buildPath, taskName + '.exe')
submitPath = os.path.join(solutionDir, hwName, taskName, 'Task', 'submit.py')

testsDir = os.path.join(solutionDir, hwName, taskName, 'Task', 'tests')

myName = 'mikhailov_nikita_m'

def runOnTest(name):
    testName = name + '.public'
    ansName = name + '.answer'
    cmd = ' '.join([solverPath, testName, ansName])
    print(cmd)
    os.system(cmd)

def submitTest(name):
    testName = name + '.public'
    testName = os.path.split(testName)[-1]
    ansName = name + '.answer'
    cmd = ' '.join(['python', submitPath, myName, taskName, testName, ansName])
    print(cmd)
    os.system(cmd)

def main():
    tests = os.listdir(testsDir)
    procs = []
    timeProcessed = time()
    for test in tests:
        name, ext = os.path.splitext(test)
        if ext != '.public':
            continue
        name = os.path.join(testsDir, name)
        print('Running ' + test)
        proc = mp.Process(target=runOnTest, args=(name,))
        procs.append(proc)
        proc.start()
    
    for proc in procs:
        proc.join()
        proc.terminate()
    timeProcessed = timeProcessed - time()
    if timeProcessed >= 10:
        print("TIME LIMIT\n")
        return
        
    procs = []
    for test in tests:
        name, ext = os.path.splitext(test)
        if ext != '.answer':
            continue
        name = os.path.join(testsDir, name)
        print('Submit ' + test)
        proc = mp.Process(target=submitTest, args=(name,))
        procs.append(proc)
        proc.start()
    
    for proc in procs:
        proc.join()
        proc.terminate()
    

if __name__ == '__main__':
    main()