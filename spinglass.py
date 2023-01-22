import numpy as np

def main():
    
    N = 1000 # sidelength of grid
    Ns = N*N # number of spins
   
    model = 1 # 0 = Ising, 1 = Antiferromagnetic Ising, 2 = Spin Glass

    numits = 100000

    h = 0
    beta = 10
    print("Preparing neighbors.")
    J = getJ(N,Ns,model)
    print("Neighbors prepared.")
    #print(np.max(J-J.T))
    #print(np.sum(np.abs(J),axis=0))

    sig = np.random.choice([-1,1],[Ns,1])

    for it in range(numits): 
        
        if it%500 == 0:
            print(f"it {it}/{numits}")
            wfile = open(f"/home/matt/Documents/mmk/ising/runs/antiferromagnet/N100/sigAF_{it}.bin","wb")
            wfile.write(sig)
            wfile.close()

        sigtrial = np.copy(sig)

        # Define the current energy state
        current_energy = getH(sig,J,h)


        spin = np.random.randint(Ns)
        sigtrial[spin]*=-1

        # Define the proposed trial energy
        trial_energy = getH(sigtrial,J,h)
        #print(current_energy, trial_energy)
        
        # Calculate the difference in energy
        delta = trial_energy - current_energy

        if delta < 0:
            sig = np.copy(sigtrial)
        else:

            # Generate a random number between 0 and 1
            random_num = np.random.rand()
            
            # Compare the random number with the acceptance probability
            if beta*delta>30:
                acceptance_prob = 0
            else:
                acceptance_prob = np.exp(-beta*delta)
            if random_num < acceptance_prob:
                # Accept the trial move
                sig = np.copy(sigtrial)



def get1Dcoord(i,j,N):
    
    return j+i*N

def get2Dcoord(k,N):
    
    i = k//N
    j = k-i*N
    return i,j

def getJ(N,Ns,model):

    J = np.zeros([Ns,Ns])

    for i in range(Ns):
        ix, iy = get2Dcoord(i,N)
        #print(f"Point = ({ix},{iy})")
        neighbors = [[ix-1,iy],[ix+1,iy],[ix,iy-1],[ix,iy+1]]
        #print(f"Before: neighbors = {neighbors}")
        for l in range(len(neighbors)):
            for n in range(len(neighbors[l])):
                if neighbors[l][n] < 0:
                    neighbors[l][n] = N-1
                elif neighbors[l][n] > N-1:
                    neighbors[l][n] = 0

        #print(f"After: neighbors = {neighbors}\n")
        neigh1D = [get1Dcoord(neigh[0],neigh[1],N) for neigh in neighbors]

        for k in neigh1D:
            J[i,k] = 1
    if model == 0:
        J *= 1
    elif model == 1:
        J *= -1
    elif model == 2:
        weights = np.random.choice([1,-1],Ns*Ns)
        k = -1
        for i in range(Ns):
            for j in range(i,Ns):
                k+=1    
                J[i,j]*=weights[k]
                J[j,i]*=weights[k]
        del weights


    return J

def getH(sig,J,h):
    return (sig.T@J@sig)[0,0]*(-0.5) + h*np.sum(sig)

if __name__=="__main__":
    main()
