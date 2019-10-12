import ctypes
import math
import subprocess
import binascii
import numpy as np
import struct
import string
import random
import os

#---------------------------------------
# Utility functions

REPO_PATH = os.path.dirname(os.path.realpath(__file__))+'/'

def writebytestringtobinfile(bytestring,filename):
	with open(filename,mode='wb') as f:
		f.write(bytestring)

def readbinfiletobytestring(filename):
	bytestring = bytearray()
	with open(filename,"rb") as l:
		bytestring = l.read()
	return bytestring

# Takes in list of integers (python int) and creates a file with these that schifra can use
def erasurelisttoerasurefile(erasurelist,filename):
	with open(filename,"wb") as g:
		for i in erasurelist:
			bytestring_error_loc = struct.pack('H', i)
			g.write(bytestring_error_loc)

def compilecppfileforRScodec(RS_code_len,RS_data_len):
	# we need to compile the cpp code with some parameters already known. hence this rigmarole.
	f = open(REPO_PATH+"/RS_paramaters_from_python.hpp", "w")
	f.write("const std::size_t code_length = 65535;\n")
	f.write("const std::size_t fec_length  =  " + str(RS_code_len-RS_data_len) +";\n")
	f.write("const std::size_t data_length = code_length - fec_length;\n")
	f.write("const std::size_t field_descriptor                =   16;\n")
	f.write("const std::size_t generator_polynomial_index      =    0;\n")
	f.write("const std::size_t generator_polynomial_root_count =  " + str(RS_code_len-RS_data_len) +";\n")
	f.close()
	subprocess.call('g++ -std=c++11 -O3 -march=native -o '+REPO_PATH+'/schifra_RS_16bit_fileio '+REPO_PATH+'/schifra_RS_16bit_fileio.cpp', shell = True)


# this function takes a bytestring of length 2*codeword_data_len (ie codeword_data_len 16bit symbols)
# and an int codeword_redundancy as input.
# It outputs an RS encoded bytestring of length 2*(codeword_data_len+codeword_redundancy)
# (ie codeword_data_len+codeword_redundancy symbols)
def RS_encode_16bit(inputbytestring,codeword_data_len,codeword_redundancy):
	RS_code_len = 65535
	RS_data_len = RS_code_len - codeword_redundancy

	# we need to compile the cpp code with some parameters already known. hence this rigmarole.
	compilecppfileforRScodec(RS_code_len,RS_data_len)
	

	# Now we pad our input bytestring.
	padding_len = RS_data_len - codeword_data_len
	RS_inputbytestring = inputbytestring.rjust(2*RS_data_len,b'0') 

	writebytestringtobinfile(RS_inputbytestring,REPO_PATH+"/trialinput.dat")
	subprocess.call([REPO_PATH+"./schifra_RS_16bit_fileio 1 "+REPO_PATH+"/trialinput.dat "+REPO_PATH+"/trialoutput.schifra 1 "+REPO_PATH+"/trialerasurelocationfile.dat"], shell = True)

	encodedbytestring = readbinfiletobytestring(REPO_PATH+"/trialoutput.schifra")

	# now puncturing the encoded codeword
	puncturedbytestring = encodedbytestring[2*padding_len:]


	# cleanup all the temporary files we have created
	os.remove(REPO_PATH+"/schifra_RS_16bit_fileio")
	os.remove(REPO_PATH+"/RS_paramaters_from_python.hpp")
	os.remove(REPO_PATH+"/trialinput.dat")
	os.remove(REPO_PATH+"/trialoutput.schifra")

	return puncturedbytestring



# this function takes a bytestring of length 2*(codeword_data_len+codeword_redundancy)
# (ie (codeword_data_len+codeword_redundancy) 16bit symbols)
# and an int codeword_redundancy as input.
# it also takes as input a list of ints of erasure locations (this can be the empty list)
# It outputs an RS decoded bytestring of length 2*(codeword_data_len)
# (ie codeword_data_len 16bit symbols)
def RS_decode_16bit(inputbytestring,codeword_data_len,codeword_redundancy,proto_erasure_loc_list):
	RS_code_len = 65535
	RS_data_len = RS_code_len - codeword_redundancy

	# we need to compile the cpp code with some parameters already known. hence this rigmarole.
	compilecppfileforRScodec(RS_code_len,RS_data_len)

	# Now we pad our input bytestring. Appending to the list of 'erasure locations' for the puncturing
	padding_len = RS_data_len - codeword_data_len
	RS_inputbytestring = inputbytestring.rjust(2*RS_code_len,b'0') 
	erasure_loc_list = []
	for x in proto_erasure_loc_list:
		erasure_loc_list.append(x+padding_len)


	writebytestringtobinfile(RS_inputbytestring,REPO_PATH+"/trialinput.dat")

	if (len(erasure_loc_list)==0):
		subprocess.call([REPO_PATH + "/schifra_RS_16bit_fileio 0 "+REPO_PATH+"/trialinput.dat "+REPO_PATH+"/trialoutput.schifra 0 "+REPO_PATH+"/trialerasurelocationfile.dat"], shell = True)
	else:
		# writing erasure_loc_list to file so schifra can use it
		erasurelisttoerasurefile(erasure_loc_list,REPO_PATH+"/trialerasurelocationfile.dat")
		subprocess.call([REPO_PATH+"/schifra_RS_16bit_fileio 0 "+REPO_PATH+"/trialinput.dat "+REPO_PATH+"/trialoutput.schifra 1 "+REPO_PATH+"/trialerasurelocationfile.dat"], shell = True)
	
	if (os.path.isfile(REPO_PATH+"/trialoutput.schifra")):
		decodedbytestring = readbinfiletobytestring(REPO_PATH+"/trialoutput.schifra")
		# now puncturing the decoded codeword
		puncturedbytestring = decodedbytestring[2*padding_len:]
		os.remove(REPO_PATH+"/trialoutput.schifra")
	else:
		puncturedbytestring = b''.rjust(2*(codeword_data_len),b'0')
	

	# cleanup all the temporary files we have created
	os.remove(REPO_PATH+"/schifra_RS_16bit_fileio")
	os.remove(REPO_PATH+"/RS_paramaters_from_python.hpp")
	os.remove(REPO_PATH+"/trialinput.dat")
	os.remove(REPO_PATH+"/trialerasurelocationfile.dat")

	return puncturedbytestring


'''
data_len = 5

length_bit_string = 16*data_len
bit_string = ''.join(random.choices(['0','1'], k = length_bit_string))
byte_string = binascii.unhexlify(((hex(int(bit_string,2)))[2:]).zfill(length_bit_string//4))


num_erasures = 4
#erasureloclist = [i for i in range(num_erasures)]
erasureloclist = [0,1,3,2]


encccc = RS_encode_16bit(byte_string,data_len,num_erasures)
zenc = encccc[:2*(num_erasures+data_len)]
tenc = zenc.ljust(2*(data_len+num_erasures),b'0')
finalbytestring = RS_decode_16bit(tenc,data_len,num_erasures,erasureloclist)

print((byte_string))
print((finalbytestring)) 
print(len(byte_string))
print(len(finalbytestring))
'''


# the input is a list of size numreads of bytestrings 
# each bytestring has length (2*symbolsperread) (ie symbolsperread symbols)
# the output is a list of size symbolsperread of bytestrings
# each bytestring has length (2*numreads) (ie numreads symbols)
def listofreadstolistofRSinputdata(listofreads):
	symbolsperread = int((len(listofreads[0]))/2)
	numreads = len(listofreads)
	
	listofRSinputdata = []
	for RSinputstring_id in range(symbolsperread):
		RSinputstring_aslist = []
		for read_id in range(numreads):
			read = listofreads[read_id]
			RSinputstring_aslist.append(read[2*RSinputstring_id:2*RSinputstring_id+2])
		RSinputstring = b''.join(RSinputstring_aslist)
		listofRSinputdata.append(RSinputstring)

	return listofRSinputdata
'''
data_len = 5

length_bit_string = 16*data_len
bit_string = ''.join(random.choices(['0','1'], k = length_bit_string))
byte_string = binascii.unhexlify(((hex(int(bit_string,2)))[2:]).zfill(length_bit_string//4))

listofreads = [byte_string,byte_string,(b'').rjust(2*data_len,b'0')]
print(listofreads)
listRS = listofreadstolistofRSinputdata(listofreads)
print(listRS)
'''

# the input is a list of size symbolsperread of bytestrings
# each bytestring has length (2*numreads) (ie numreads symbols)
# the output is a list of size numreads of bytestrings 
# each bytestring has length (2*symbolsperread) (ie symbolsperread symbols)
def listofRSoutputdatatolistofreads(listofRSoutputdata):
	numreads = int((len(listofRSoutputdata[0]))/2)
	symbolsperread = len(listofRSoutputdata)
	
	listofreads = []
	for read_id in range(numreads):
		readstring_aslist = []
		for RSoutputstring_id in range(symbolsperread):
			RSinputstring = listofRSoutputdata[RSoutputstring_id]
			readstring_aslist.append(RSinputstring[2*read_id:2*read_id+2])
		readstring = b''.join(readstring_aslist)
		listofreads.append(readstring)

	return listofreads

'''
data_len = 5

length_bit_string = 16*data_len
bit_string = ''.join(random.choices(['0','1'], k = length_bit_string))
byte_string = binascii.unhexlify(((hex(int(bit_string,2)))[2:]).zfill(length_bit_string//4))

listofreads = [byte_string,byte_string,(b'').rjust(2*data_len,b'0')]
print(listofreads)
listRS = listofreadstolistofRSinputdata(listofreads)
print(listRS)
outputreads = listofRSoutputdatatolistofreads(listRS)
print(outputreads)
'''






# this function takes in list of [readid,corruptedread]
# corrupted read is a bytestring of length (2*symbolsperread) (ie symbolsperread symbols)
# we are also given as input int totalnumreads which is the number of 
# reads in uor list that we would get if there were no erasures.
# size of input list is \leq this integer
# output is twofold
# first output is a list of bytestrings which are RS decoder input data
# each bytestring has length (2*totalnumreads) (ie totalnumreads symbols)
# the second output is a list of integers containing erasure locations
# (in the erased location of bytestring, we insert dummy strings)
def listofcorruptedreadstoRSinputdataanderasureloclist(listofcorruptedreads,totalnumreads):
	symbolsperread = int((len(listofcorruptedreads[0][1]))/2)
	dummyread = b''.rjust(2*symbolsperread,b'0')
	listofcorruptedreadswithdummies = [dummyread for i in range(totalnumreads)]
	erasure_loc_list = [i for i in range(totalnumreads)]
	for j in listofcorruptedreads:
		listofcorruptedreadswithdummies[j[0]] = j[1]
		erasure_loc_list.remove(j[0])
	
	listofRSinputdata = listofreadstolistofRSinputdata(listofcorruptedreadswithdummies)
	
	return listofRSinputdata,erasure_loc_list

'''
print("testing stuff")	
listofcorruptedreads = [[2,byte_string],[0,byte_string]]
print(listofcorruptedreads)
listofRSinputdata,erasure_loc_list = listofcorruptedreadstoRSinputdataanderasureloclist(listofcorruptedreads,4)

print(erasure_loc_list)
print(listofRSinputdata)
'''
# End of utility functions
#---------------------------------------


# the input is a list of size numreads of bytestrings 
# each bytestring has length (2*symbolsperread) (ie symbolsperread symbols)
# we also input int redundancy which is the amount of redundancy per RS codeword
# the output is a list of size (numreads+redundancy) of RS encoded (vertically) bytestrings 
# each bytestring has length (2*symbolsperread) (ie symbolsperread symbols)
def MainEncoder(listofreads,redundancy):
	numreads = len(listofreads)
	
	listofRSinputdata = listofreadstolistofRSinputdata(listofreads)

	listofRSoutputdata = []
	for RSinputdata in listofRSinputdata:
		RSoutputdata = RS_encode_16bit(RSinputdata,numreads,redundancy)
		listofRSoutputdata.append(RSoutputdata)

	listofRSencodedreads = listofRSoutputdatatolistofreads(listofRSoutputdata)
	return listofRSencodedreads


# this function takes in list of [readid,corruptedread]
# corrupted read is a bytestring of length (2*symbolsperread) (ie symbolsperread symbols)
# we also input int redundancy which is the amount of redundancy per RS codeword
# we get as input totalnumreads which is the number of
# reads in our list that we would get if there were no erasures.
# size of input list is \leq this integer
# we also input int redundancy which is the amount of redundancy per RS codeword
# the input is a list of size (numreads+redundancy) of RS encoded (vertically) bytestrings 
# each bytestring has length (2*symbolsperread) (ie symbolsperread symbols)
def MainDecoder(listofcorruptedreads,redundancy,totalnumreads):
	listofRSinputdata,erasure_loc_list = listofcorruptedreadstoRSinputdataanderasureloclist(listofcorruptedreads,totalnumreads)
	
	outputnumreads = totalnumreads - redundancy
	listofRSoutputdata = []
	for RSinputdata in listofRSinputdata:
		RSoutputdata = RS_decode_16bit(RSinputdata,outputnumreads,redundancy,erasure_loc_list)
		listofRSoutputdata.append(RSoutputdata)

	listofRSdecodedreads = listofRSoutputdatatolistofreads(listofRSoutputdata)
	return listofRSdecodedreads

'''
data_len = 8

length_bit_string = 16*data_len
bit_string = ''.join(random.choices(['0','1'], k = length_bit_string))
byte_string = binascii.unhexlify(((hex(int(bit_string,2)))[2:]).zfill(length_bit_string//4))

redundancy = 6
listofreads = [byte_string,byte_string,(b'').rjust(2*data_len,b'0')]
totalnumreads = len(listofreads)


enclist = MainEncoder(listofreads,redundancy)

totalnumreads = len(enclist)


enclistwithid = []
numreadshere = len(enclist)
for i in range(numreadshere):
	enclistwithid.append([i,enclist[i]])

# introducing erasures
numerasures = redundancy - 4
enclistwithid = enclistwithid[:-numerasures]

#introducing errors
numerrors = 5
for i in range(numerrors):
	enclistwithid[i][1] = b''.ljust(len(enclistwithid[0][1]),b'0')

declist = MainDecoder(enclistwithid,redundancy,totalnumreads)
print(listofreads)
print(declist)
'''








