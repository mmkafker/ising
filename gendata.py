import numpy as np

N = 1000
Ns = N*N
numtrials = 1000000000

pfile = open("flips.bin","wb")
flips = np.random.randint(Ns,size=numtrials,dtype=np.int32)
pfile.write(flips)
#print(flips)
pfile.close()


probfile = open("probs.bin","wb")
probs = np.random.rand(numtrials)
probfile.write(probs)
#print(probs)
probfile.close()


spinfile = open("spins.bin","wb")
spins = np.random.choice(np.array([1, -1],dtype=np.int32), Ns)
spinfile.write(spins)
#print(spins)
spinfile.close()


