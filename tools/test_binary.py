import struct
from numpy.random import rand

f=open('/Users/dmelgar/misc/test_binary.blurg','wb')
# Write an integer
f.write(struct.pack('i',362))
# Write a float (single precision 4byte,32bits)
f.write(struct.pack('f',-3.123e4))
# Write an array of doubles (double precision, 8bytes,64bits)
d=rand(10)*1000  #Some random numbers
f.write(struct.pack('d'*len(d),*d))
# Write shitty string
string='poop'
f.write(struct.pack('s'*len(string),*string))

#Donezo, close
f.close()

f=open('/Users/dmelgar/misc/test_binary.blurg','rb')
bs=f.read()