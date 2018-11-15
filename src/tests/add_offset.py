import numpy as np

a1 = np.ones((4,5))
a2 = np.zeros((4,5))
a2[0,0]=1
a2[:,3]=2

offset = (2,1)

print(a1)
print(a2)
print(offset)

print('a1+a2')
print(a1+a2)

nullar = np.zeros((len(a1)+offset[0], len(a1[0])+offset[1]))
print(nullar)
nullar[:-offset[0], :-offset[1]] += a1
a1 = np.copy(nullar)
print(a1)

a1[offset[0]:, offset[1]:] += a2
print(a1)
