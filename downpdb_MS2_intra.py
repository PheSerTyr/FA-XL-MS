from selenium import webdriver
from selenium.webdriver.common.by import By
import requests
from Bio import SeqIO
import re
from Bio.PDB import PDBParser

distancecutoff=20

def compute_ca_distance(pdb_file):
    structure = PDBParser().get_structure('pdb', pdb_file)
    model = structure[0]
    k_residues = []
    r_residues = []
    chain=model["A"]
    residue_shift = -2
    for residue in chain:
        if residue_shift==-2:
            residue_shift = residue.get_id()[1]-1
        if residue.get_resname() == 'LYS':
            k_residues.append(residue)
        elif residue.get_resname() == 'ARG':
            r_residues.append(residue)
    distances = {}
    for k_residue in k_residues:
        for r_residue in r_residues:
            try:
                k_ca = k_residue['CA'].get_coord()
                r_ca = r_residue['CA'].get_coord()
                distance = ((k_ca[0]-r_ca[0])**2+(k_ca[1]-r_ca[1])**2+(k_ca[2]-r_ca[2])**2)**0.5
                distances[(k_residue.get_id()[1]-residue_shift,r_residue.get_id()[1]-residue_shift)]=[distance,False]
                #print(distance)
                '''print(f'C alpha distance between K{str(k_residue.get_id()[1])} '
                      f'and R{str(r_residue.get_id()[1])}: {distance} angstroms')'''
            except KeyError:
                pass
    return distances

def retrieve_c_alpha(pdb_file, chain_id, residue_number):
    # Create a PDB parser object
    parser = PDBParser()

    # Parse the PDB file
    structure = parser.get_structure("structure", pdb_file)

    # Get the model from the structure (assuming single model)
    model = structure[0]

    # Get the specified chain
    chain = model[chain_id]

    for residue in chain:
        residue_shift = residue.get_id()[1]-1
        break

    # Get the residue at the specified position
    residue = chain[residue_number+residue_shift]

    b_factors = [atom.get_bfactor() for atom in residue]
    confidence = sum(b_factors) / len(b_factors)

    # Get the C-alpha atom
    c_alpha = residue['CA']

    # Retrieve the coordinates of the C-alpha atom
    c_alpha_point = c_alpha.get_coord()

    return [c_alpha_point,confidence]

def retrieve_string_in_parentheses(input_string):
    pattern = r'\((.*?)\)'  # Regex pattern to match strings within parentheses
    matches = re.findall(pattern, input_string)  # Find all matches in the input string
    
    return matches

def find_peptide_positions(protein_sequence, peptide_sequence):
    positions = []
    peptide_length = len(peptide_sequence)
    protein_length = len(protein_sequence)

    # Iterate over the protein sequence
    for i in range(protein_length - peptide_length + 1):
        # Check if the current substring matches the peptide sequence
        if protein_sequence[i:i + peptide_length] == peptide_sequence:
            positions.append(i)

    return positions


def retrieve_protein_sequence(pdb_file):
    sequences = []
    with open(pdb_file, "r") as handle:
        for record in SeqIO.parse(handle, "pdb-atom"):
            if record.seq:
                sequences.append(str(record.seq))
    return sequences

def download_file(url, save_path):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(save_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)

        print("File downloaded successfully!")
        
    except requests.exceptions.RequestException as e:
        print("Error occurred while downloading the file:", e)


path="D:\\MS_DATA\\plus_20250211\\Brain_intra_trypsin_rep3.tsv"



pdbpath="D:\\MS_DATA\\AlphaPdb\\"


f=open(path,"r")
sw=open(path+"_distance.tsv","w")

#"P02769":"E:\\FA_BSA\\pdb4f5s.ent"
pdb={"P02769":"D:\\MS_DATA\\distance\\4f5s.pdb"} #若有需要不下载的pdb文件在这里指定"P02769 uniprot序号":"pdb文件路径"

distances=compute_ca_distance("D:\\MS_DATA\\distance\\4f5s.pdb")

sw.write(f.readline().strip("\n")+"\n")

for line in f.readlines():
    words=line.split("\t")
    #try:
        #if float(words[43])>0.1 or float(words[45])>0.1: #前后两个cutoff
            #continue
    #except:
        #continue

    uniprotids=list(set(words[9].split("|")) & set(words[22].split("|")))
    pep1=words[11]
    link1=int(retrieve_string_in_parentheses(words[12])[-1])
    pep2=words[24]
    link2=int(retrieve_string_in_parentheses(words[25])[-1])
    d=[]

    nofileflag=True
    for uniprotid in uniprotids:
        print(uniprotid)
        if uniprotid not in pdb.keys():
            # Set up Selenium with Chrome WebDriver
            driver = webdriver.Chrome()

            url = "https://alphafold.ebi.ac.uk/entry/"+uniprotid

            # Navigate to the URL
            driver.get(url)

            try:
                # Find the PDB file button element by its XPath
                pdb_button = driver.find_element(By.XPATH, '//a[contains(text(), "PDB file")]')

                # Retrieve the href attribute value (link) of the PDB file button
                pdb_url = pdb_button.get_attribute('href')
                print("PDB File Link:", pdb_url)
                save_path = pdbpath + pdb_url.split("/")[-1]
                try:
                    download_file(pdb_url, save_path)
                    pdb[uniprotid]=save_path
                except:
                    print("download failed")
                    continue
                # Close the browser
                driver.quit()
            except:
                print("PDB file link not found.")
                # Close the browser
                driver.quit()
                continue
        else:
            print("already downloaded")

        # Provide the path to your PDB file
        pdb_file = pdb[uniprotid]

        sequences = retrieve_protein_sequence(pdb_file)
        for seq in sequences:
            print("Protein Sequence:")
            print(seq)
        
        # Provide the protein sequence and peptide sequence
        protein_sequence = sequences[0]

        try:
            position1 = find_peptide_positions(protein_sequence, pep1)[0]
            position2 = find_peptide_positions(protein_sequence, pep2)[0]
        except:
            sw.write(line.strip("\n")+"\t")
            sw.write("N/A")
            continue
        
        s=[position1+link1,position2+link2]
        site1=s[0]
        site2=s[1]

        p1,p1c=retrieve_c_alpha(pdb_file, "A", s[0])
        p2,p2c=retrieve_c_alpha(pdb_file, "A", s[1])

        print(p1)
        print(p2)
        print(((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)**0.5)

        if tuple(s) in distances.keys():
            distances[tuple(s)][1]=True
        elif (s[1],s[0]) in distances.keys():
            distances[(s[1],s[0])][1]=True


        #if p1c>70 and p2c>70:
        if p1c>0 and p2c>0:
            d.append(((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)**0.5)
        else:
            nofileflag=False
    sw.write(line.strip("\n")+"\t")
    if len(d)!=0:
        sw.write("{0:.2f}".format(min(d)))
    elif nofileflag:
        sw.write("N/A")
    else:
        sw.write("low confidence")
    sw.write("\n")

sw.close()

for distancecutoff in [5,10,15,20,25,30,40,50]:
    print(distancecutoff)

    tst=[0,0,0,0] #TP,FP,TN,FN

    for k,v in distances.items():
        if v[0]<distancecutoff and v[1]:
            tst[0]+=1
        elif v[0] <distancecutoff:
            tst[3]+=1
        elif v[1]:
            tst[1]+=1
        else:
            tst[2]+=1


    print(tst)


print("end")
