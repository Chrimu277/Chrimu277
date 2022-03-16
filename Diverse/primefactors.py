# This is a sample Python script.

import pandas as pd
import matplotlib as mpl
import numpy as np

# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

_smallp = ( 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
    313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
    409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
    499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
    601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
    691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
    907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997)

def primefactors(n):

    factors = []
    n_dash = n

#check for primefactors and dividability

    while n_dash != 1:
        for i in range(0,len(_smallp)):

            if n_dash%_smallp[i] == 0:
                n_dash = n_dash/_smallp[i];
                factors.append(_smallp[i])
                break;

    return factors



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Strg+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print("This program gives you prime factors of given numbers in a file")

    pf_file = open("primefactors_file.txt")
    pf_file_text = pf_file.readlines();     #list of text
    nr_lines = len(pf_file_text)

    numbers_array = np.array([0 for i in range(nr_lines)],dtype=int) # make array as long as
    numbers_list = []

    for i in range(0,nr_lines):
        numbers_array[i] = int(pf_file_text[i])
        numbers_list.append(int(pf_file_text[i]))

    #calculate primefactors
    for i in range(0, nr_lines):
        print("Prime factors of "+str(numbers_array[i])+" are: "+str(primefactors(numbers_array[i])))

    

