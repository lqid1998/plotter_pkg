def file_read(path):
    with open(path) as file_object: 
        contents = file_object.readlines()
    return contents

def get_charge_arr(contents):
    charge_arr = []
    line_temp = []
    contents_len=len(contents)
    for i in range(contents_len):
        if i > 1 and i < (contents_len-4):
            line_temp = contents[i].split()
            charge_arr.append(float(line_temp[4]))
    return charge_arr

def get_charge_sum(path, start, end):
    contents = file_read(path)
    charge_arr = get_charge_arr(contents)
    charge_sum = 0
    for i in range(start, end):
        atom_charge = charge_arr[i]
        charge_sum = charge_sum + charge_arr[i]
    return charge_sum
