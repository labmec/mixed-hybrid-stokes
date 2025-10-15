import json


viscosity = 0.0005
Re = int(1/viscosity)
size = [1, 2, 3, 4, 5, 6]

for k in size:
    # json Files
    File_json = "/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen.json"
    newFile_json = f"TestOseen_Re{Re}_k{k}.json"

    with open(File_json, "r") as f:
        data = json.load(f)

    data["MeshName"] = f"/home/itopo/NS_Research/mixed-hybrid-stokes/TestOseen_k{k}"
    data["Domain"][0]["viscosity"] = viscosity

    with open(newFile_json, "w") as f:
        json.dump(data, f, indent=4)



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
