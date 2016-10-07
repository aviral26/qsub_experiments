import subprocess

l = []
for i in (0, 1):
    l.append(subprocess.Popen(['python', 'DummySubprocess.py'], shell=False, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE))

print "Waiting for {0} processes.".format(len(l))

outerrs = [p.communicate("This is a random line.") for p in l]

print "Outerrs are " + str(outerrs)
