import subprocess
string="AAAAAAAA BBBBBBB CCCCCC"
var=string.split(" ")
print(var)
print(var[0])
output = subprocess.check_output(
        [ r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\muscle5.1.win64.exe",
        "-align", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Trial.txt",
        "-output", r"E:\\SeniorI_Fall2022_CUFE_HEM\\Bioinformatics\\Assignments\\Multiple Sequence Alignment\\MyTask\\Trial__output.fasta"],
        text=True)

# ofile = open("my_fasta.txt", "w")

# for i in range(len(var)):

# ofile.write(">" + list_name[i] + "\n" +list_seq[i] + "\n")

# #do not forget to close it

# ofile.close()
#======================================

