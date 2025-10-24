import json


viscosity = 0.000333333333
Re = int(1/viscosity)
size = [1, 2, 3, 4, 5, 6]

# for k in size:
#     # json Files
#     #File_json = "/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen.json"
#     File_json = "/home/marina/programming/NS-Research/mixed-hybrid-stokes/TestOseen.json"
#     newFile_json = f"TestOseen_Re{Re}_k{k}.json"

#     with open(File_json, "r") as f:
#         data = json.load(f)

#     #data["MeshName"] = f"/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen_k{k}"
#     data["MeshName"] = f"/home/marina/programming/NS-Research/mixed-hybrid-stokes/TestOseen_k{k}"
#     data["Domain"][0]["viscosity"] = viscosity

#     with open(newFile_json, "w") as f:
#         json.dump(data, f, indent=4)



# geo Files
# File_geo = "/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen.geo"
# newFile_geo = f"/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen_k{k}.geo"

# with open(File_geo, "r", encoding="utf-8") as f:
#     lines = f.readlines()

# nElements = int(2**k)
# lines[3] = f"el = {nElements};\n"

# with open(newFile_geo, "w", encoding="utf-8") as f:
#     f.writelines(lines)

# # close all files
# f.close()


# Writing errors in .dat format
List_errors = []
for k in size:
    # File_txt = f"/home/itopo/NS_Research/mixed-hybrid-stokes/Results_TestOseen/TestOseen_Re{Re}_k{k}_Errors.txt"
    File_txt = f"/home/marina/programming/NS-Research/mixed-hybrid-stokes/Results_TestOseen/TestOseen_Re{Re}_k{k}_Errors.txt"

    with open(File_txt, "r", encoding="utf-8") as f:
        lines = f.readlines()

    List_errors.append(lines[2])

List_errors_dat = []
List_errors_mathematica = []
List_errors_latex = []
new_List_errors_latex = []
index = 0
for line in List_errors:
    if 'e' in line:
        line = line.replace('e', 'E')
        List_errors_dat.append(line)
for line in List_errors:
    line = line.replace(' ', ', ')
    line = line.replace('e', '*10^')
    List_errors_mathematica.append(line)
for line in List_errors:
    List_errors_latex.append(line.split())
for line in List_errors_latex:
    aux = []
    for i, error in enumerate(line):
        if i in [0, 2, 4, 6]:
            aux.append('$'+error+'$')
    new_List_errors_latex.append(aux)


with open(f"Results_TestOseen_Re{Re}.txt", "w", encoding="utf-8") as f:
    for line in List_errors_dat:
        f.write(line)
        f.write('\n')
    f.write('\n')
    for line in List_errors_mathematica:
        f.write(line)
        f.write('\n')
    f.write('\n')
    for line in new_List_errors_latex:
        f.write(' & '.join(line))
        f.write('\n')

# print(List_errors_latex)
# print(new_List_errors_latex)