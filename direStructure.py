import os
import shutil

def doit():
    shutil.rmtree('c:/pythontest', ignore_errors=True)
    os.mkdir('./pythontest')
    os.chdir('./pythontest')

    for i in range(0,3):
        fileName = 'folder%d' % i
        print fileName
        os.mkdir(fileName)
        os.chdir(fileName)
        for j in range(0,3):
            fileName = 'folder%d_%d' % (i,j)
            print fileName
            os.mkdir(fileName)
            os.chdir(fileName)
            for k in range(0,3):
                try:
                    with open('file%d_%d_%d.txt' % (i,j,k), 'w'):
                        pass
                except IOError:
                    pass
            os.chdir('..')
        os.chdir('..')

