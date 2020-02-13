import struct

# The input file MUST meet the following criteria:
  # 1) All numbers must be on a single line
  # 2) All numbers must be separated by only a sinlge white space (\s)
  # 3) The number of data points is the first number; N
  # 4) The number of attributes/dimensions, for each point, is the second number; M
  # 5) The numbers that follow must be the amount of N*M

inputFile = "./particles_0_64_ascii" # change this to the file you want to convert
# inputFile = "./part64.bin" # change this to the file you want to convert
outPutFile = './part64.bin' # Change this to wht you want the output to be named

with open(inputFile, mode='rb') as file:
  # fileContent = file.read()
  # print((len(struct.unpack("i" * (len(fileContent) // 4), fileContent)) - 2)/3)
  # Read the entire file, and split it into a list, based off of the white space
  fileContent = file.read().split(' ')

  with open(outPutFile, mode='wb') as f:
    # The first 2 numbers have to be converted as integers, they are read as ints by the code
    f.write(struct.pack("i" * 2, *[int(i) for i in fileContent[:2]]))
    # The rest of the numbers are converted to floats, as they are read in the code
    f.write(struct.pack("f" * (len(fileContent) - 2), *[float(i) for i in fileContent[2:]]))
    
# Exit cleanly
exit(0)
