# This is a set of scripts to edit files
# These functions were written with colvars in mind

# clean file will remove all lines which contain a particular string
def clean_file(file_name,remove_words):

    import subprocess

    for word in remove_words:
        subprocess.call(["sed", "-i",  '/'+word+'/d', file_name])


# round_numbers will go through all the numbers in a file and replace them with a given precision (or number of decimal places)
def round_numbers(input_file,output_name,n):

    fout = open(output_name, 'w')

    with open(input_file, 'r') as fin:
        for line in fin:
            values = line.split()
            for val in values:
                new_val = round_to_n(float(val), n)  # n digit precision
                fout.write(str(new_val))
                fout.write(' ')
            fout.write('\n')

# round a number to n significant figures
def round_to_n(x,n):
    from math import log10, floor
    from numpy import sign
    if x==0:
        new_val=0
    else:
        new_val = round(x, -int(floor(log10(abs(x)))) + (n - 1))
    return new_val


# will take a file and break each column in to a separate file
def break_into_columns(input_file,column_labels):

    with open(input_file, 'r') as fin:
        x=fin.readline()
        str_list = list(filter(None, x.split(' ')))
        str_list = [x for x in str_list if x != '\n']
        columns=len(str_list)
        for i in range(columns):
            with open(input_file, 'r') as fin:
                tmp=open(column_labels+'_'+str(i), 'w')
                for line in fin:
                    xtmp=(list(filter(None, line.split(' '))))
                    tmp.write(xtmp[i].rstrip('\r\n'))
                    tmp.write('\n')
