

irc_forward = open("irc_forward.xyz", "r") #opens forward direction of irc file from psi4

irc_backward = open("irc_backward.xyz", "r") #opens forward direction of irc file from psi4

count = 0

natoms = int(irc_forward.readline())

block_size = 2+natoms
total_number_of_lines = sum(1 for line in irc_forward) + 1
irc_forward.close()
num_structs_f = int(total_number_of_lines/block_size)
total_number_of_lines = sum(1 for line in irc_backward) + 1
irc_backward.close()
num_structs_b = int(total_number_of_lines/block_size)
irc_forward = open("irc_forward.xyz", "r") #opens forward direction of irc file from psi4
forward = irc_forward.readlines()
irc_backward = open("irc_backward.xyz", "r") #opens forward direction of irc file from psi4
backward = irc_backward.readlines()

count = 0
new_irc_for = open("irc_forward_new.xyz", "w")
for i in range(int(round(num_structs_f/5))):
    for j in range(block_size):
        new_irc_for.write(forward[j+count])
    count += 5*block_size
new_irc_for.close()

count = 0
new_irc_back = open("irc_backward_new.xyz", "w")
for i in range(int(round(num_structs_b/5))):
    for j in range(block_size):
        new_irc_back.write(backward[j+count])
    count += 5*block_size
new_irc_back.close()

count = 0
new_irc_for = open("irc_forward_new.xyz", "r")
new_irc_strip = open("irc_forward_new_strip.xyz", "w") 
for line in new_irc_for:
    if not line.strip():
        new_irc_strip.write("IRC point %d\n" %(count))
        count = count+1 
    else:
        new_irc_strip.write(line)
new_irc_for.close()
new_irc_strip.close()
count = 0

new_irc_back = open("irc_backward_new.xyz", "r")
new_irc_strip = open("irc_backward_new_strip.xyz", "w")
for line in new_irc_back:
    if not line.strip():
        new_irc_strip.write("IRC point %d\n" %(-1*count))
        count = count+1
    else:
        new_irc_strip.write(line)
new_irc_back.close()
new_irc_strip.close()
