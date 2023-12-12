import filecmp
 
f1 = r"C:\Users\20212774\OneDrive - TU Eindhoven\Desktop\Boeken\Y3\Q2\OGO comp bio\GitHub\storingDataInDNA\test_data\Check.txt"
f2 = r"C:\Users\20212774\OneDrive - TU Eindhoven\Desktop\Boeken\Y3\Q2\OGO comp bio\GitHub\storingDataInDNA\test_data\output.txt"
 
# shallow comparison
result = filecmp.cmp(f1, f2)
print(result)
# deep comparison
result = filecmp.cmp(f1, f2, shallow=False)
print(result)

file1 = open(f1,'r')
file2 = open(f2,'r')

file1_lines = file1.readlines()
file2_lines = file2.readlines()

for i in range(len(file1_lines)):
    if file1_lines[i] != file2_lines[i]:
        print("Line " + str(i+1) + " doesn't match.")
        print("------------------------")
        print("File1: " + file1_lines[i])
        print("File2: " + file2_lines[i])

file1.close()
file2.close()