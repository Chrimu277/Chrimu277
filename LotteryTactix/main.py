import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

#6 aus 45
N_pick = 6
N_total = 45

def test():
    # mein lotto pick
    mypick = randomPickWithRules(N_pick,N_total,checkRules = checkRules)
    print("my pick is \n",mypick)

    # lottery lotto pick
    lotterypick = randomLotteryPick(N_pick,N_total)
    print("lottery pick is \n", lotterypick)

    # check equals
    N_equals,equals = getEquals(mypick,lotterypick)
    print(equals,"(",N_equals,") equal")

def getEquals(list1,list2):
    n_equals = 0
    equals = []
    for i in range(len(list1)):
        for j in range(len(list2)):
            if list1[i] == list2[j]:
                n_equals += 1
                equals.append(list1[i])
                break;
    return n_equals,equals

def checkRules(picks):
    return True;

def randomLotteryPick(N_p,N_t):

    arr = []
    while len(arr) < N_p:
        rand = np.random.randint(1,N_t+1,size=1)
        if not arr.__contains__(rand):
            arr.append((rand[0]))

    return arr

def randomPickWithRules(N_p,N_t,checkRules:staticmethod):
    valid = False

    arr = []
    while valid is False:
        arr.clear()
        while len(arr) < N_p:
            rand = np.random.randint(1, N_t + 1, size=1)
            if not arr.__contains__(rand):
                arr.append((rand[0]))
        valid = checkRules(arr)
    return arr

if __name__ == '__main__':
    np.random.seed(dt.now().microsecond)
    tstart = dt.now()

    N_draws = 10000 # 1000 games ~~ 20 jahre lotto
    N_tipps = 1
    equals_matrix = np.zeros((N_pick + 1, N_draws)) # spiele x tipps x [0,1,2,3,4,5,6]

    # several games
    for i in range(N_draws):
        lotterypick = randomLotteryPick(N_pick,N_total)

        # several draws per game
        for j in range(N_tipps):
            mypick = randomPickWithRules(N_pick,N_total,checkRules=checkRules)
            n_equals, equals = getEquals(mypick, lotterypick)

            if(n_equals == 6):
                print("hit a 6 -- ")
                print(i,j)
            print(mypick)


            equals_matrix[n_equals,i] += 1

    for e in range(N_pick+1):
        plt.plot(np.arange(N_draws), equals_matrix[e])
    plt.legend(["0","1","2","3","4","5","6"])
    plt.show()

    # print results
    for i in range(N_pick+1):
        percentage:float = 100 * np.sum(equals_matrix[i]) / (N_draws * N_tipps)
        print("total", i,"equals found [%]:", percentage)
        if percentage > 0:
            print(" ( 1 in ", int(np.round(100/percentage)))
        else:
            print("no chance..")
