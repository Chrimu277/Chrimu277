# Author Mühleder Christoph
# Date April 4 2021
# DONE

import numpy as np

#--------ZIP arrays-----------------------------------------

names = ['Peter Parker','Clark Kent','Wide Wilson']
heroes = ['Peter Parker','Superman','Deadpool']
universe = ['Marvel1','DC','Marvel2']

for na, he, un in zip(names,heroes,universe):
    print(f'{na} is actually {he} from the {un} Universe')


#--- numpy vector/array  -------------

list1 = [10,20,30,40,50]
print("len(list1) is ", len(list1))
print("list1 is ", list(list1))

# --- strings to numbers in List  -----
str_nums = ["4", "8", "6", "5", "3", "2", "8", "9", "2", "5"]
int_nums = map(int, str_nums)
print("int num array is", list(int_nums))

# --- list to array ----

array2 = np.asarray(list1)
print("array2 is ",array2)
print("len(array2) is ", len(array2))

# --- make arrays np
array3 = np.arange(1,100,2)
print("len(array3) is ",len(array3))
print("array3 is ",array3)

array4 = np.linspace(start=10, stop = 50, num= 2000, endpoint=True, dtype=float)
print("beginning of array4: ", array4[:-1]) # [:-1] first to 2nd last element
print("end of array4: ",array4[0:25]) # [0:25] first to 25th  element

array5 = array4[::100] # [::100] takes every 100th position
array6 = array5[::-1]   # reverse list
print("len(array5) is ",len(array5))
print("a5 vs a6 : ", array5, "\n and \n", array6)


