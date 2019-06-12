import multiprocessing as mp
import os

solutionDir = 'C:\\DOPT-2019-Practice'
hwName = 'hw2'
taskName = 'tsp'
buildPath = os.path.join(solutionDir, hwName,'x64\\Release\\')

solverPath = os.path.join(buildPath, taskName + '.exe')

testsDir = os.path.join(solutionDir, hwName, taskName, 'Task', 'tests')

def runOnTest(name):
    testName = name + '.public'
    ansName = name + '.answer'
    cmd = solverPath + ' ' + testName + ' ' + ansName
    os.system(cmd)



def main():
    tests = os.listdir(testsDir)
    procs = []
    for test in tests:
        name, ext = os.path.splitext(test)
        if ext != '.public':
            continue

        print('Running ' + test)
        proc = mp.Process(target=runOnTest, args=(name,))
        procs.append(proc)
        proc.start()
    
    for proc in procs:
        proc.join()
        

if __name__ == '__main__':
    main()