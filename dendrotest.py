import os
import dendropy
import argparse
from Bio import SeqIO

#Читаем параметры командной строки
descriptionText = "Produces phylogenetic trees annotated with aminoacid substitutions"
parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-fasta", dest="aln", required=True, help="Path to alignment FASTA file")
parser.add_argument("-engine", dest="engine", required=False, help="Phylogenetic tree building engine (raxml, raxml-ng or iqtree)", default="raxml")
parser.add_argument("-raxml", dest="raxml", required=False, help="Path to RaXML executable", default="raxmlHPC")
parser.add_argument("-fasttree", dest="fasttree", required=False, help="Path to FastTree executable", default="FastTreeMP")
parser.add_argument("-raxml-ng", dest="raxmlng", required=False, help="Path to RaXML-NG executable", default="raxml-ng")
parser.add_argument("-iqtree", dest="iqtree", required=False, help="Path to IQtree executable", default="iqtree")
parser.add_argument("-orfs", dest="orfs", required=False, help="Path to file with open reading frames' coordinates", default="0")
args = parser.parse_args()

#Пути к исполняемым файлам
RAXML_PATH = args.raxml
FASTTREE_PATH = args.fasttree
RAXMLNG_PATH = args.raxmlng
IQTREE_PATH = args.iqtree
if not os.path.isfile(args.raxml):
    print("Warning! RaXML executable not found")
if not os.path.isfile(args.fasttree):
    print("Warning! FastTree executable not found")
if not os.path.isfile(args.raxmlng):
    print("Warning! RaXML-NG executable not found")
if not os.path.isfile(args.iqtree):
    print("Warning! IQtree executable not found")

#Читаем входное выравнивание, выбрасываем дублирующиеся по названию последовательности
seqs = []
seq_names = {}
print("Open input sequence file")
for seq_record in SeqIO.parse(args.aln, 'fasta'):
	tmp_str = seq_record.id.upper()
	if tmp_str not in seq_names.keys():
		seq_names[tmp_str] = 1
		seqs.append(seq_record)
print("Write reduced sequence file")
with open("aln.fas", 'w') as w:
	for seq in seqs:
		w.write(">"+seq.id+"\n"+str(seq.seq)+"\n")
print('Total used: '+str(len(seqs))+' sequences')

#Строим дерево, выбираем метод в зависимости от числа последовательностей
if len(seqs) <= 500:
    if args.engine == 'raxml':
        print("We have less than 500 sequences.\nStarting RAxML tree search using GTRGAMMA model and 1000 fast bootstrap replications\n\n")
        os.system(RAXML_PATH+" -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -n topol -T "+str(os.cpu_count())+" -s aln.fas -o "+seqs[0].id)
        tree_name = "RAxML_bipartitions.topol"
        print("RAxML tree search finished\n\n---------------------------------------------------\n\n")
elif len(seqs) > 500 and len(seqs) <= 1000:
	print("We have more than 500 sequences and less than 1000 sequences.\nStarting RAxML tree search using GTRGAMMA model and 10 tree search attempts without bootstrap\n\n")
	os.system(RAXML_PATH+" -m GTRGAMMA -p 12345 -# 10 -n topol -T "+str(os.cpu_count())+" -s aln.fas -o "+seqs[0].id)
	tree_name = "RAxML_bestTree.topol"
	print("RAxML tree search finished\n\n---------------------------------------------------\n\n")
elif len(seqs) > 1000 and len(seqs) <= 3000:
	print("We have more than 1000 sequences and less than 3000 sequences.\nStarting RAxML tree search using GTRCAT model and 1 tree search attempt without bootstrap\n\n")
	os.system(RAXML_PATH+" -m GTRCAT -p 12345 -n topol -T "+str(os.cpu_count())+" -s aln.fas -o "+seqs[0].id)
	tree_name = "RAxML_bestTree.topol"
	print("RAxML tree search finished\n\n---------------------------------------------------\n\n")
elif len(seqs) > 3000 and len(seqs) <= 10000:
	print("We have more than 3000 sequences and less than 10000 sequences.\nStarting combined FastTree and RAxML tree search using GTRCAT model and 1 tree search attempt without bootstrap\n\n")
	if not os.path.isfile("FastTree_unrooted.topol"):
		os.system(FASTTREE_PATH+" -nt -gtr -gamma < aln.fas > FastTree_unrooted.topol")
	tmp_tree = dendropy.Tree.get(path="FastTree_unrooted.topol", schema="newick")
	for node in tmp_tree.leaf_node_iter():
		tmp = str(node.taxon)
		tmpp = tmp.replace(" ","_")
		node.taxon.label = tmpp[1:len(tmpp)-1]
	#root_node = tmp_tree.find_node_with_taxon_label("'"+seqs[0].id+"'")
	#tmp_tree.to_outgroup_position(root_node)
	tmp_tree.ladderize(ascending=False)
	tmp_tree.resolve_polytomies()
	tmp_tree.write(path="FastTree_unrooted.topol", schema="newick")
	os.system(RAXML_PATH+" -m GTRCAT -n topol -T "+str(os.cpu_count())+" -s aln.fas -t FastTree_unrooted.topol -o "+seqs[0].id)
	tree_name = "RAxML_bestTree.topol"
	print("RAxML tree search finished\n\n---------------------------------------------------\n\n")
elif len(seqs) > 10000:
	print("We have more than 10000 sequences.\nStarting FastTree tree search \n\n")
	if not os.path.isfile("FastTree_unrooted.topol"):
		os.system(FASTTREE_PATH+" -nt -gtr -gamma < aln.fas > FastTree_unrooted.topol")
	tmp_tree = dendropy.Tree.get(path="FastTree_unrooted.topol", schema="newick")
	for node in tmp_tree.leaf_node_iter():
		tmp = str(node.taxon)
		tmpp = tmp.replace(" ","_")
		node.taxon.label = tmpp[1:len(tmpp)-1]
		#print(node.label)
	root_node = tmp_tree.find_node_with_taxon_label("'"+seqs[0].id+"'")
	#print(str(seqs[0].id))
	#print(root_node.label)
	tmp_tree.to_outgroup_position(root_node)
	#tmp = tmp_tree.as_string(schema="newick").split(" ")
	tmp_tree.ladderize(ascending=False)
	tmp_tree.resolve_polytomies()
	#tmp = tmp_tree.as_string(schema="newick").replace("'","")
	#with open("FastTree_rooted.topol", 'w') as w:
	#	w.write(tmp)
	tmp_tree.write(path="FastTree_rooted.topol", schema="newick")
	tree_name = "FastTree_rooted.topol"
	print("FastTree tree search finished\n\n---------------------------------------------------\n\n")

#Реконструкция последовательностей
print("Starting RAxML ancestral state reconstruction\n\n")
os.system("raxmlHPC -m GTRGAMMA -f A -n ances -s aln.fas -T "+str(os.cpu_count())+" -t "+tree_name)
#./raxml-ng-static --ancestral --msa tst37H3seq.fas --tree fasttree37_sdr.tree --model GTR+G --prefix tst8
print("RAxML ancestral state reconstruction finished\n\n---------------------------------------------------\n\n")

#Открываем оба дерева
print("Opening tree files")
tree = dendropy.Tree.get(path=tree_name, schema="newick")
anctree = dendropy.Tree.get(path="RAxML_nodeLabelledRootedTree.ances", schema="newick")

#Упорядочиваем оба дерева
print("Sorting trees")
tree.ladderize(ascending=False)
anctree.ladderize(ascending=False)

#Переносим бутстреп в аннотации
cou = 0
print("\n")
for node in tree.postorder_node_iter():
	if node.label != None:
		node.annotations['bootstrap'] = node.label
		node.label = None
		print("Moving bootstrap value for node "+str(cou)+" of "+str(len(seqs)-1)+"\r", end="")
		cou +=1

#Переносим названия листьев в аннотации
cou = 0
print("\n")
for node in tree.leaf_node_iter():
	#print(node.taxon)
	tax_tmp = str(node.taxon)
	tax = tax_tmp.replace(" ","_")
	node.annotations['name'] = tax[1:len(tax)-1]
	node.label = str(node.taxon)
	print("Moving names to annotations for leaf "+str(cou)+" of "+str(len(seqs))+"\r", end="")
	cou +=1

#Тест вывода деревьев в консоль
#tree.print_plot()
#anctree.print_plot()


##1.Первый проход - копирование номеров внутренних узлов с второго на первое дерево от листьев к корню. 

cou = 0
print("\n")
for nodeanc in anctree.postorder_node_iter():
	if nodeanc.label == None and nodeanc.taxon != None:
		anc_tax = str(nodeanc.taxon)
		#print(anc_tax[1:len(anc_tax)-1])
		nodeancparent = nodeanc.parent_node
		node = tree.find_node_with_label(str(nodeanc.taxon))
		#print(node.taxon)
		nodeparent = node.parent_node
		nodeparent.annotations['node_number'] = nodeancparent.label
		nodeparent.label = nodeancparent.label
		#print(nd.parent_node)
		#print(node)
		#print(ndp.label)
	else:
		#print(nodeanc.label)
		nodeancparent = nodeanc.parent_node
		node = tree.find_node_with_label(str(nodeanc.label))
		if node == None:
			child_anc = nodeanc.child_nodes()
			child_anc_first = child_anc[0]
			node_child = tree.find_node_with_label(str(child_anc_first.label))
			node = node_child.parent_node
			node.label = nodeanc.label
			node.annotations['node_number'] = node.label
		#print(node.label)
		nodeparent = node.parent_node
		if nodeparent != None:
			nodeparent.annotations['node_number'] = nodeancparent.label
	print("Moving name to label for node "+str(cou)+"\r", end="")
	cou +=1

##2.Второй проход - нанесение замен на ветви и листья
#Открываем файл с реконструированными последовательностями
print("\n")
print("Open reconstructed sequences file")
with open('RAxML_marginalAncestralStates.ances', 'r') as f:
	anc_seq_arr = f.readlines()
with open('RAxML_marginalAncestralStates.ances.fas', 'w') as w:
	for seq in anc_seq_arr:
		seq_arr = seq.split(' ')
		w.write(">"+seq_arr[0]+"\n"+seq_arr[1]+"\n")
for seq_record in SeqIO.parse('RAxML_marginalAncestralStates.ances.fas', 'fasta'):
	seqs.append(seq_record)
print('Total: '+str(len(seqs))+' sequences')

#Сохраняем все нуклеотидные последовательности в файл
print("Saving all nucleotide sequences")
with open('nucl_seqs_tst.fas', 'w') as w:
	for el in seqs:
		w.write(">"+el.id+"\n"+str(el.seq)+"\n")

#Транслируем последовательности
prot_seq = {}
cou = 0
print("\n")
for seq in seqs:
	seq_str_tmp = str(seq.seq)
	seq_str = seq_str_tmp.upper()
	trans_seq_str = ""
	i = 0
	while i < len(seq_str)-1:
		#print(seq_str[i:i+3])
		if seq_str[i:i+3] in ["ATG"]:
			trans_seq_str += "M"
		elif seq_str[i:i+3] in ["GCT", "GCC", "GCA", "GCG", "GCN", "GCR", "GCY", "GCM", "GCK", "GCS", "GCW", "GCH", "GCB", "GCV", "GCD"]:
			trans_seq_str += "A"
		elif seq_str[i:i+3] in ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "CGN", "CGR", "CGY", "CGM", "CGK", "CGS", "CGW", "CGH", "CGB", "CGV", "CGD", "AGR"]:
			trans_seq_str += "R"
		elif seq_str[i:i+3] in ["AAT", "AAC", "ATY"]:
			trans_seq_str += "N"
		elif seq_str[i:i+3] in ["GAT", "GAC", "GAY"]:
			trans_seq_str += "D"
		elif seq_str[i:i+3] in ["TGT", "TGC", "TGY"]:
			trans_seq_str += "C"
		elif seq_str[i:i+3] in ["CAA", "CAG", "CAR"]:
			trans_seq_str += "Q"
		elif seq_str[i:i+3] in ["GAA", "GAG", "GAR"]:
			trans_seq_str += "E"
		elif seq_str[i:i+3] in ["GGT", "GGC", "GGA", "GGG", "GGN", "GCR", "GGY", "GGM", "GGK", "GGS", "GGW", "GGH", "GGB", "GGV", "GGD"]:
			trans_seq_str += "G"
		elif seq_str[i:i+3] in ["CAT", "CAC", "CAY"]:
			trans_seq_str += "H"
		elif seq_str[i:i+3] in ["ATT", "ATC", "ATA"]:
			trans_seq_str += "I"
		elif seq_str[i:i+3] in ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "TTR", "CTN", "CTR", "CTY", "CTM", "CTK", "CTS", "CTW", "CTH", "CTB", "CTV", "CTD"]:
			trans_seq_str += "L"
		elif seq_str[i:i+3] in ["AAA", "AAG", "AAR"]:
			trans_seq_str += "K"
		elif seq_str[i:i+3] in ["TTT", "TTC", "TTY"]:
			trans_seq_str += "F"
		elif seq_str[i:i+3] in ["CCT", "CCC", "CCA", "CCG", "CCN", "CCR", "CCY", "CCM", "CCK", "CCS", "CCW", "CCH", "CCB", "CCV", "CCD"]:
			trans_seq_str += "P"
		elif seq_str[i:i+3] in ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "AGY", "TCN", "TCR", "TCY", "TCM", "TCK", "TCS", "TCW", "TCH", "TCB", "TCV", "TCD"]:
			trans_seq_str += "S"
		elif seq_str[i:i+3] in ["ACT", "ACC", "ACA", "ACG", "ACN", "ACR", "ACY", "ACM", "ACK", "ACS", "ACW", "ACH", "ACB", "ACV", "ACD"]:
			trans_seq_str += "T"
		elif seq_str[i:i+3] in ["TGG"]:
			trans_seq_str += "W"
		elif seq_str[i:i+3] in ["TAT", "TAC", "TAY"]:
			trans_seq_str += "Y"
		elif seq_str[i:i+3] in ["GTT", "GTC", "GTA", "GTG", "GTN", "GTR", "GTY", "GTM", "GTK", "GTS", "GTW", "GTH", "GTB", "GTV", "GTD"]:
			trans_seq_str += "V"
		elif seq_str[i:i+3] in ["TAG", "TGA", "TAA", "TRA", "TAR"]:
			trans_seq_str += "*"
		elif seq_str[i:i+3] in ["---", "???"]:
			trans_seq_str += "-"
		else:
			trans_seq_str += "?"
		i +=3
	prot_seq[seq.id] = trans_seq_str
	print("Translating sequence "+str(cou)+" of "+str(len(seqs))+"\r", end="")
	cou +=1
#print(prot_seq.keys())

#Сохраняем транслированные последовательности в файл
print("\n")
print("Saving translated sequences")
with open('prot_seqs_tst.fas', 'w') as w:
	for el in prot_seq.keys():
		w.write(">"+el+"\n"+prot_seq[el]+"\n")


#Наносим замены
cou = 0
print("\n")
for node in tree.postorder_node_iter():
	if node.parent_node != None:
		nodeparent = node.parent_node
		tax_tmp = str(node.label)
		tax = tax_tmp.replace(" ","_")
		#print(tax)
		tax_parent_tmp = str(nodeparent.label)
		tax_parent = tax_parent_tmp.replace(" ","_")
		#print(tax_parent)
		if node.is_leaf():
			seq2 = prot_seq[tax[1:len(tax)-1]]
			seq1 = prot_seq[tax_parent]
			subst_str = ""
			for i in range(len(seq1)):
				if seq1[i] != seq2[i]:
					subst = seq1[i]+str(i+1)+seq2[i]
					if len(subst_str) != 0:
						subst_str += ", "
					subst_str +=subst
			if len(subst_str) != 0:
				node.annotations['name_and_subs'] = '"'+tax[1:len(tax)-1]+" ["+subst_str+"]"+'"'
			else:
				node.annotations['name_and_subs'] = '"'+tax[1:len(tax)-1]+'"'
		else:
			if node.child_nodes() != None:
				children = node.child_nodes()
				tax_child1_tmp = str(children[0].label)
				tax_child1 = tax_child1_tmp.replace(" ","_")
				tax_child2_tmp = str(children[1].label)
				tax_child2 = tax_child2_tmp.replace(" ","_")
				seq2 = prot_seq[tax]
				seq1 = prot_seq[tax_parent]
				if children[0].is_leaf():
					seq3 = prot_seq[tax_child1[1:len(tax_child1)-1]]
				else:
					seq3 = prot_seq[tax_child1]
				if children[1].is_leaf():
					seq4 = prot_seq[tax_child2[1:len(tax_child2)-1]]
				else:
					seq4 = prot_seq[tax_child2]
				subst_str = ""
				for i in range(len(seq1)):
					if seq1[i] != seq2[i] and seq2[i] == seq3[i] and seq2[i] == seq4[i]:
						subst = seq1[i]+str(i+1)+seq2[i]
						if len(subst_str) != 0:
							subst_str += ", "
						subst_str +=subst
				edg = node.edge
				if len(subst_str) != 0:
					edg.annotations['subs'] = '"'+"["+subst_str+"]"+'"'
	print("Placing substitutions on branch "+str(cou)+" of "+str(len(seqs)*2)+"\r", end="")
	cou +=1


#Тест сохранения в NEXUS для Figtree
print("\n")
print("Saving result tree")
tree.write(path="RAxML_final.nex", schema="nexus")
