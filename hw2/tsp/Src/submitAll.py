import multiprocessing as mp
import os

solutionDir = 'C:\\DOPT-2019-Practice'
hwName = 'hw2'
taskName = 'tsp'
buildPath = os.path.join(solutionDir, hwName,'x64\\Release\\')

solverPath = os.path.join(buildPath, taskName + '.exe')
submitPath = os.path.join(solutionDir, hwName, taskName, 'Task', 'submit.py')

testsDir = os.path.join(solutionDir, hwName, taskName, 'Task', 'tests')

myName = 'test_name'

def runOnTest(name):
    testName = name + '.public'
    ansName = name + '.answer'
    cmd = ' '.join([solverPath, testName, ansName])
    print(cmd)
    os.system(cmd)

def submitTest(name):
    testName = name + '.public'
    ansName = name + '.answer'
    cmd = ' '.join(['python', submitPath, myName, taskName, testName, ansName])
    print(cmd)
    os.system(cmd)

def main():
    tests = os.listdir(testsDir)
    procs = []
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
    

if __name__ == '__main__':
    main()