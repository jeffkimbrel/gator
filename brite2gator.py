import json

kegg = json.loads(open("/Users/kimbrel1/Desktop/working/amino_acid_transport/ko02000_filtered.json", "r").read())

pathways = {}

brite = kegg['name']
for level1 in kegg['children']:
   high = level1['name']
   for level2 in level1['children']:
      mid = level2['name']
      for level3 in level2['children']:
         low = level3['name']
         pathways[low] = {}
         # print(f"--- {high}; {mid}; {low} ---")
         if 'children' in level3:
            for level4 in level3['children']:
               line = level4['name']
               KO = line.split().pop(0)
               # print(f"{KO}***{line}")
               pathways[low][KO] = line

for path, data in pathways.items():

   paths = []

   for KO, description in data.items():
      gene_name,description = description.split("; ")
      gene_name = gene_name.replace(", ", "/")
      gene_name = gene_name.replace("  ", " ")
      gene_name = gene_name.replace(" ", "_")

      paths.append(gene_name)

      # for genes
      # print(f"{gene_name}\t{KO}\t{description}\t{path}")

   # for paths
   print(f"{path}\t{' -> '.join(paths)}")
