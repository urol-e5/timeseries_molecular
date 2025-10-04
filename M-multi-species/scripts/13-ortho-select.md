# 13-ortho-select


``` bash
head ../output/11-orthology-analysis/ortholog_groups.csv
```

    group_id,apul,peve,ptua,type,avg_identity
    OG_00001,FUN_000185-T1,Peve_00037402,Pocillopora_meandrina_HIv1___RNAseq.g28886.t1,three_way,71.059
    OG_00002,FUN_000189-T1,Peve_00038462,Pocillopora_meandrina_HIv1___RNAseq.g28888.t3,three_way,68.72583333333334
    OG_00003,FUN_000190-T1,Peve_00038463,Pocillopora_meandrina_HIv1___RNAseq.g28889.t1,three_way,55.1975
    OG_00004,FUN_000191-T1,Peve_00038464,Pocillopora_meandrina_HIv1___RNAseq.g28890.t1,three_way,72.047
    OG_00005,FUN_000192-T1,Peve_00038466,Pocillopora_meandrina_HIv1___RNAseq.g28893.t1,three_way,81.15866666666666
    OG_00006,FUN_000193-T1,Peve_00038467,Pocillopora_meandrina_HIv1___RNAseq.g28894.t1,three_way,70.21933333333334
    OG_00007,FUN_000194-T1,Peve_00038470,Pocillopora_meandrina_HIv1___RNAseq.g28898.t1,three_way,77.11833333333334
    OG_00008,FUN_000195-T1,Peve_00038471,Pocillopora_meandrina_HIv1___RNAseq.g28899.t1,three_way,68.38
    OG_00009,FUN_000197-T1,Peve_00027182,Pocillopora_meandrina_HIv1___RNAseq.g28901.t1,three_way,56.498333333333335

``` bash
wc -l ../output/11-orthology-analysis/ortholog_groups.csv
```

    18327 ../output/11-orthology-analysis/ortholog_groups.csv

## Annotated

``` bash
head ../output/12-ortho-annot/ortholog_groups_annotated.csv
```

    group_id,apul,peve,ptua,type,avg_identity,query,accession,id,reviewed,protein_name,organism,pident,length,evalue,bitscore,title,go_ids,go_bp,go_cc,go_mf,goslim_ids,goslim_names
    OG_00001,FUN_000185-T1,Peve_00037402,Pocillopora_meandrina_HIv1___RNAseq.g28886.t1,three_way,71.059,FUN_000185-T1,Q32LQ0,AMPE_BOVIN,,Glutamyl aminopeptidase (EAP) (EC 3.4.11.7) (Aminopeptidase A) (AP-A) (CD antigen CD249),Bos taurus (Bovine),41.5,881.0,1.4e-197,691.4,sp|Q32LQ0|AMPE_BOVIN Glutamyl aminopeptidase OS=Bos taurus OX=9913 GN=ENPEP PE=2 SV=1,GO:0003081; GO:0004177; GO:0004230; GO:0005615; GO:0005737; GO:0005886; GO:0006508; GO:0008217; GO:0008270; GO:0008283; GO:0016477; GO:0042277; GO:0043171; GO:0070006,cell migration [GO:0016477]; cell population proliferation [GO:0008283]; peptide catabolic process [GO:0043171]; proteolysis [GO:0006508]; regulation of blood pressure [GO:0008217]; regulation of systemic arterial blood pressure by renin-angiotensin [GO:0003081],cytoplasm [GO:0005737]; extracellular space [GO:0005615]; plasma membrane [GO:0005886],aminopeptidase activity [GO:0004177]; glutamyl aminopeptidase activity [GO:0004230]; metalloaminopeptidase activity [GO:0070006]; peptide binding [GO:0042277]; zinc ion binding [GO:0008270],GO:0003824; GO:0005615; GO:0005886; GO:0016787; GO:0048870; GO:0050886; GO:0140096,"catalytic activity; extracellular space; plasma membrane; hydrolase activity; cell motility; endocrine process; catalytic activity, acting on a protein"
    OG_00002,FUN_000189-T1,Peve_00038462,Pocillopora_meandrina_HIv1___RNAseq.g28888.t3,three_way,68.72583333333334,FUN_000189-T1,Q9DCT6,,,,,38.0,171.0,4.6e-18,92.8,sp|Q9DCT6|BAP18_MOUSE BPTF-associated chromatin complex component 1 OS=Mus musculus OX=10090 GN=Bacc1 PE=1 SV=1,,,,,,
    OG_00003,FUN_000190-T1,Peve_00038463,Pocillopora_meandrina_HIv1___RNAseq.g28889.t1,three_way,55.1975,,,,,,,,,,,,,,,,,
    OG_00004,FUN_000191-T1,Peve_00038464,Pocillopora_meandrina_HIv1___RNAseq.g28890.t1,three_way,72.047,FUN_000191-T1,Q95108,,,,,44.1,136.0,3.3e-29,129.4,"sp|Q95108|THIOM_BOVIN Thioredoxin, mitochondrial OS=Bos taurus OX=9913 GN=TXN2 PE=1 SV=2",,,,,,
    OG_00005,FUN_000192-T1,Peve_00038466,Pocillopora_meandrina_HIv1___RNAseq.g28893.t1,three_way,81.15866666666666,FUN_000192-T1,Q1LWX3,,,,,54.3,188.0,3.5e-53,209.5,sp|Q1LWX3|NTAQ1_DANRE Protein N-terminal glutamine amidohydrolase OS=Danio rerio OX=7955 GN=ntaq1 PE=2 SV=1,,,,,,
    OG_00006,FUN_000193-T1,Peve_00038467,Pocillopora_meandrina_HIv1___RNAseq.g28894.t1,three_way,70.21933333333334,FUN_000193-T1,Q5I0K7,,,,,53.2,158.0,3.4e-41,169.5,sp|Q5I0K7|ALG13_RAT UDP-N-acetylglucosamine transferase subunit ALG13 OS=Rattus norvegicus OX=10116 GN=Alg13 PE=1 SV=1,,,,,,
    OG_00007,FUN_000194-T1,Peve_00038470,Pocillopora_meandrina_HIv1___RNAseq.g28898.t1,three_way,77.11833333333334,FUN_000194-T1,P70182,,,,,60.2,447.0,3.4e-147,523.5,sp|P70182|PI51A_MOUSE Phosphatidylinositol 4-phosphate 5-kinase type-1 alpha OS=Mus musculus OX=10090 GN=Pip5k1a PE=1 SV=2,,,,,,
    OG_00008,FUN_000195-T1,Peve_00038471,Pocillopora_meandrina_HIv1___RNAseq.g28899.t1,three_way,68.38,FUN_000195-T1,B8D9S0,,,,,38.5,148.0,9e-19,95.1,sp|B8D9S0|FMT_BUCA5 Methionyl-tRNA formyltransferase OS=Buchnera aphidicola subsp. Acyrthosiphon pisum (strain 5A) OX=563178 GN=fmt PE=3 SV=1,,,,,,
    OG_00009,FUN_000197-T1,Peve_00027182,Pocillopora_meandrina_HIv1___RNAseq.g28901.t1,three_way,56.498333333333335,FUN_000197-T1,Q8TCT7,,,,,44.7,488.0,1.6e-89,331.6,sp|Q8TCT7|SPP2B_HUMAN Signal peptide peptidase-like 2B OS=Homo sapiens OX=9606 GN=SPPL2B PE=1 SV=2,,,,,,

``` bash
wc -l ../output/12-ortho-annot/ortholog_groups_annotated.csv
```

    18373 ../output/12-ortho-annot/ortholog_groups_annotated.csv

# Physiology based selection

## Filter for specific GO terms

``` bash
# Filter orthologs for metabolic and stress-related processes
echo "Filtering for:"
echo "- Glycolysis (GO:0006096)"
echo "- Gluconeogenesis (GO:0006094)"
echo "- Lipolysis/lipid catabolism (GO:0016042)"
echo "- Fatty acid beta oxidation (GO:0006635)"
echo "- Starvation (GO:0042594)"
echo "- Lipid biosynthesis (GO:0008610)"
echo "- Protein catabolic process (GO:0030163)"
echo "- Reproductive Process (GO:0022414)"
echo ""

# Create output directory if it doesn't exist
mkdir -p ../output/13-ortho-select

# Extract header
head -1 ../output/12-ortho-annot/ortholog_groups_annotated.csv > ../output/13-ortho-select/selected_orthologs.csv

# Filter for the specified GO terms and reproduction-related processes
grep -E 'GO:0006096|GO:0006094|GO:0016042|GO:0006635|GO:0042594|GO:0008610|GO:0030163|GO:0022414' \
  ../output/12-ortho-annot/ortholog_groups_annotated.csv \
  >> ../output/13-ortho-select/selected_orthologs.csv

echo "Total orthologs selected:"
wc -l ../output/13-ortho-select/selected_orthologs.csv
```

    Filtering for:
    - Glycolysis (GO:0006096)
    - Gluconeogenesis (GO:0006094)
    - Lipolysis/lipid catabolism (GO:0016042)
    - Fatty acid beta oxidation (GO:0006635)
    - Starvation (GO:0042594)
    - Lipid biosynthesis (GO:0008610)
    - Protein catabolic process (GO:0030163)
    - Reproductive Process (GO:0022414)

    Total orthologs selected:
    787 ../output/13-ortho-select/selected_orthologs.csv

``` bash
# Show first few selected orthologs
head ../output/13-ortho-select/selected_orthologs.csv
```

    group_id,apul,peve,ptua,type,avg_identity,query,accession,id,reviewed,protein_name,organism,pident,length,evalue,bitscore,title,go_ids,go_bp,go_cc,go_mf,goslim_ids,goslim_names
    OG_00044,FUN_000272-T1,Peve_00027798,Pocillopora_meandrina_HIv1___RNAseq.g28985.t1,three_way,70.59733333333334,FUN_000272-T1,O75031,HSF2B_HUMAN,,Heat shock factor 2-binding protein,Homo sapiens (Human),43.2,190.0,5.9e-32,139.0,sp|O75031|HSF2B_HUMAN Heat shock factor 2-binding protein OS=Homo sapiens OX=9606 GN=HSF2BP PE=1 SV=1,GO:0005654; GO:0005694; GO:0005829; GO:0006366; GO:0007141; GO:0007144; GO:0007283; GO:1990918,double-strand break repair involved in meiotic recombination [GO:1990918]; female meiosis I [GO:0007144]; male meiosis I [GO:0007141]; spermatogenesis [GO:0007283]; transcription by RNA polymerase II [GO:0006366],chromosome [GO:0005694]; cytosol [GO:0005829]; nucleoplasm [GO:0005654],,GO:0005654; GO:0005694; GO:0005829; GO:0006281; GO:0006351; GO:0022414; GO:0043226; GO:0140013,nucleoplasm; chromosome; cytosol; DNA repair; DNA-templated transcription; reproductive process; organelle; meiotic nuclear division
    OG_00086,FUN_000437-T1,Peve_00011193,Pocillopora_meandrina_HIv1___RNAseq.g2323.t1,three_way,77.881,FUN_000437-T1,F1MNN4,FBXW7_BOVIN,,F-box/WD repeat-containing protein 7 (F-box and WD-40 domain-containing protein 7),Bos taurus (Bovine),36.4,195.0,1.9999999999999998e-26,122.5,sp|F1MNN4|FBXW7_BOVIN F-box/WD repeat-containing protein 7 OS=Bos taurus OX=9913 GN=FBXW7 PE=1 SV=3,GO:0005634; GO:0005654; GO:0005694; GO:0005737; GO:0006281; GO:0006974; GO:0007062; GO:0010629; GO:0010992; GO:0019005; GO:0030332; GO:0031146; GO:0031398; GO:0031625; GO:0042752; GO:0042802; GO:0043130; GO:0043161; GO:0050816; GO:0050821; GO:0070374; GO:0070534; GO:0097027; GO:1901524; GO:1901800; GO:1903378; GO:1903955; GO:1990452; GO:1990756; GO:2000060; GO:2001205,DNA damage response [GO:0006974]; DNA repair [GO:0006281]; negative regulation of gene expression [GO:0010629]; negative regulation of osteoclast development [GO:2001205]; positive regulation of ERK1 and ERK2 cascade [GO:0070374]; positive regulation of oxidative stress-induced neuron intrinsic apoptotic signaling pathway [GO:1903378]; positive regulation of proteasomal protein catabolic process [GO:1901800]; positive regulation of protein targeting to mitochondrion [GO:1903955]; positive regulation of protein ubiquitination [GO:0031398]; positive regulation of ubiquitin-dependent protein catabolic process [GO:2000060]; proteasome-mediated ubiquitin-dependent protein catabolic process [GO:0043161]; protein K63-linked ubiquitination [GO:0070534]; protein stabilization [GO:0050821]; regulation of circadian rhythm [GO:0042752]; regulation of mitophagy [GO:1901524]; SCF-dependent proteasomal ubiquitin-dependent protein catabolic process [GO:0031146]; sister chromatid cohesion [GO:0007062]; ubiquitin recycling [GO:0010992],chromosome [GO:0005694]; cytoplasm [GO:0005737]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]; Parkin-FBXW7-Cul1 ubiquitin ligase complex [GO:1990452]; SCF ubiquitin ligase complex [GO:0019005],cyclin binding [GO:0030332]; identical protein binding [GO:0042802]; phosphothreonine residue binding [GO:0050816]; ubiquitin binding [GO:0043130]; ubiquitin protein ligase binding [GO:0031625]; ubiquitin-like ligase-substrate adaptor activity [GO:1990756]; ubiquitin-protein transferase activator activity [GO:0097027],GO:0005634; GO:0005654; GO:0005694; GO:0006281; GO:0030163; GO:0043226; GO:0060090; GO:0098772,nucleus; nucleoplasm; chromosome; DNA repair; protein catabolic process; organelle; molecular adaptor activity; molecular function regulator activity
    OG_00107,FUN_000499-T1,Peve_00024155,Pocillopora_meandrina_HIv1___RNAseq.g2414.t1,three_way,82.06683333333332,FUN_000499-T1,Q14432,PDE3A_HUMAN,,"cGMP-inhibited 3',5'-cyclic phosphodiesterase 3A (EC 3.1.4.17) (Cyclic GMP-inhibited phosphodiesterase A) (CGI-PDE A) (cGMP-inhibited cAMP phosphodiesterase) (cGI-PDE)",Homo sapiens (Human),37.6,1034.0,2.2e-143,511.5,"sp|Q14432|PDE3A_HUMAN cGMP-inhibited 3',5'-cyclic phosphodiesterase 3A OS=Homo sapiens OX=9606 GN=PDE3A PE=1 SV=3",GO:0001556; GO:0004114; GO:0004115; GO:0004119; GO:0005829; GO:0006629; GO:0007186; GO:0009410; GO:0016020; GO:0019933; GO:0019934; GO:0030284; GO:0040020; GO:0043066; GO:0043116; GO:0043117; GO:0046872; GO:0047555; GO:0060282; GO:0060700; GO:0071321; GO:0071560; GO:0097190; GO:0099130; GO:0106072; GO:0141162,apoptotic signaling pathway [GO:0097190]; cAMP-mediated signaling [GO:0019933]; cellular response to cGMP [GO:0071321]; cellular response to transforming growth factor beta stimulus [GO:0071560]; cGMP-mediated signaling [GO:0019934]; G protein-coupled receptor signaling pathway [GO:0007186]; lipid metabolic process [GO:0006629]; negative regulation of adenylate cyclase-activating G protein-coupled receptor signaling pathway [GO:0106072]; negative regulation of apoptotic process [GO:0043066]; negative regulation of cAMP/PKA signal transduction [GO:0141162]; negative regulation of vascular permeability [GO:0043116]; oocyte maturation [GO:0001556]; positive regulation of oocyte development [GO:0060282]; positive regulation of vascular permeability [GO:0043117]; regulation of meiotic nuclear division [GO:0040020]; regulation of ribonuclease activity [GO:0060700]; response to xenobiotic stimulus [GO:0009410],cytosol [GO:0005829]; membrane [GO:0016020],"3',5'-cGMP-inhibited cyclic-nucleotide phosphodiesterase activity [GO:0004119]; 3',5'-cyclic-AMP phosphodiesterase activity [GO:0004115]; 3',5'-cyclic-GMP phosphodiesterase activity [GO:0047555]; 3',5'-cyclic-nucleotide phosphodiesterase activity [GO:0004114]; estrogen binding [GO:0099130]; metal ion binding [GO:0046872]; nuclear estrogen receptor activity [GO:0030284]",GO:0003013; GO:0003824; GO:0005829; GO:0006629; GO:0016787; GO:0022414; GO:0060089; GO:0140110,circulatory system process; catalytic activity; cytosol; lipid metabolic process; hydrolase activity; reproductive process; molecular transducer activity; transcription regulator activity
    OG_00142,FUN_000582-T1,Peve_00034702,Pocillopora_meandrina_HIv1___TS.g6501.t1,three_way,59.576166666666666,FUN_000582-T1,D2X8K2,PA2_CONGI,,Phospholipase A2 A2-actitoxin-Cgg2a (A2-AITX-Cgg2a) (EC 3.1.1.4) (CgPLA2) (Phosphatidylcholine 2-acylhydrolase),Condylactis gigantea (Giant Caribbean anemone) (Condylactis passiflora),43.7,119.0,3.5e-25,116.3,sp|D2X8K2|PA2_CONGI Phospholipase A2 A2-actitoxin-Cgg2a OS=Condylactis gigantea OX=47073 PE=1 SV=1,GO:0005509; GO:0005543; GO:0005576; GO:0006644; GO:0016042; GO:0042151; GO:0044398; GO:0044470; GO:0044522; GO:0047498; GO:0050482; GO:0090729,arachidonate secretion [GO:0050482]; lipid catabolic process [GO:0016042]; phospholipid metabolic process [GO:0006644]; venom-mediated edema in another organism [GO:0044398]; venom-mediated myocyte killing in another organism [GO:0044522]; venom-mediated suppression of blood coagulation [GO:0044470],extracellular region [GO:0005576]; nematocyst [GO:0042151],calcium ion binding [GO:0005509]; calcium-dependent phospholipase A2 activity [GO:0047498]; phospholipid binding [GO:0005543]; toxin activity [GO:0090729],GO:0003824; GO:0005576; GO:0006629; GO:0008289; GO:0016787; GO:0043226; GO:0090729,catalytic activity; extracellular region; lipid metabolic process; lipid binding; hydrolase activity; organelle; toxin activity
    OG_00192,FUN_000898-T1,Peve_00009741,Pocillopora_meandrina_HIv1___RNAseq.g29449.t1,three_way,80.12233333333333,FUN_000898-T1,O43903,GAS2_HUMAN,,Growth arrest-specific protein 2 (GAS-2),Homo sapiens (Human),42.6,256.0,1.3e-48,194.9,sp|O43903|GAS2_HUMAN Growth arrest-specific protein 2 OS=Homo sapiens OX=9606 GN=GAS2 PE=1 SV=1,GO:0001544; GO:0001547; GO:0001725; GO:0005829; GO:0005884; GO:0006915; GO:0008017; GO:0008360; GO:0008593; GO:0016020; GO:0030728; GO:0051015; GO:0051726; GO:0051764; GO:0071711,actin crosslink formation [GO:0051764]; antral ovarian follicle growth [GO:0001547]; apoptotic process [GO:0006915]; basement membrane organization [GO:0071711]; initiation of primordial ovarian follicle growth [GO:0001544]; ovulation [GO:0030728]; regulation of cell cycle [GO:0051726]; regulation of cell shape [GO:0008360]; regulation of Notch signaling pathway [GO:0008593],actin filament [GO:0005884]; cytosol [GO:0005829]; membrane [GO:0016020]; stress fiber [GO:0001725],actin filament binding [GO:0051015]; microtubule binding [GO:0008017],GO:0005829; GO:0008092; GO:0012501; GO:0022414; GO:0030198,cytosol; cytoskeletal protein binding; programmed cell death; reproductive process; extracellular matrix organization
    OG_00195,FUN_000905-T1,Peve_00009733,Pocillopora_meandrina_HIv1___RNAseq.g29456.t1,three_way,66.77466666666666,FUN_000905-T1,P70627,FOLH1_RAT,,Glutamate carboxypeptidase 2 (EC 3.4.17.21) (Folate hydrolase 1) (Folylpoly-gamma-glutamate carboxypeptidase) (FGCP) (Glutamate carboxypeptidase II) (GCPII) (Membrane glutamate carboxypeptidase) (mGCP) (N-acetylated-alpha-linked acidic dipeptidase I) (NAALADase I) (Prostate-specific membrane antigen homolog) (Pteroylpoly-gamma-glutamate carboxypeptidase),Rattus norvegicus (Rat),37.3,753.0,5.9e-130,466.5,sp|P70627|FOLH1_RAT Glutamate carboxypeptidase 2 OS=Rattus norvegicus OX=10116 GN=Folh1 PE=1 SV=1,GO:0004180; GO:0004181; GO:0005886; GO:0006508; GO:0006760; GO:0008233; GO:0009986; GO:0016020; GO:0016805; GO:0030163; GO:0043065; GO:0046872; GO:1904492; GO:1904493,folic acid-containing compound metabolic process [GO:0006760]; positive regulation of apoptotic process [GO:0043065]; protein catabolic process [GO:0030163]; proteolysis [GO:0006508],cell surface [GO:0009986]; membrane [GO:0016020]; plasma membrane [GO:0005886],Ac-Asp-Glu binding [GO:1904492]; carboxypeptidase activity [GO:0004180]; dipeptidase activity [GO:0016805]; metal ion binding [GO:0046872]; metallocarboxypeptidase activity [GO:0004181]; peptidase activity [GO:0008233]; tetrahydrofolyl-poly(glutamate) polymer binding [GO:1904493],GO:0003824; GO:0005886; GO:0006575; GO:0016787; GO:0030163; GO:0140096,"catalytic activity; plasma membrane; modified amino acid metabolic process; hydrolase activity; protein catabolic process; catalytic activity, acting on a protein"
    OG_00207,FUN_000926-T1,Peve_00009718,Pocillopora_meandrina_HIv1___RNAseq.g29467.t1,three_way,63.42333333333334,FUN_000926-T1,Q9R0N8,SYT6_MOUSE,,Synaptotagmin-6 (Synaptotagmin VI) (SytVI),Mus musculus (Mouse),29.5,319.0,9.1e-33,142.9,sp|Q9R0N8|SYT6_MOUSE Synaptotagmin-6 OS=Mus musculus OX=10090 GN=Syt6 PE=1 SV=2,GO:0001786; GO:0005544; GO:0005829; GO:0005886; GO:0007340; GO:0009898; GO:0016020; GO:0030672; GO:0042802; GO:0042803; GO:0046872; GO:0046982; GO:0048306; GO:0097038,acrosome reaction [GO:0007340],cytoplasmic side of plasma membrane [GO:0009898]; cytosol [GO:0005829]; membrane [GO:0016020]; perinuclear endoplasmic reticulum [GO:0097038]; plasma membrane [GO:0005886]; synaptic vesicle membrane [GO:0030672],calcium-dependent phospholipid binding [GO:0005544]; calcium-dependent protein binding [GO:0048306]; identical protein binding [GO:0042802]; metal ion binding [GO:0046872]; phosphatidylserine binding [GO:0001786]; protein heterodimerization activity [GO:0046982]; protein homodimerization activity [GO:0042803],GO:0005829; GO:0005886; GO:0008289; GO:0022414,cytosol; plasma membrane; lipid binding; reproductive process
    OG_00237,FUN_000992-T1,Peve_00004455,Pocillopora_meandrina_HIv1___RNAseq.g29194.t1,three_way,92.068,FUN_000992-T1,P82471,GNAQ_RAT,,Guanine nucleotide-binding protein G(q) subunit alpha (EC 3.6.5.-) (Guanine nucleotide-binding protein alpha-q),Rattus norvegicus (Rat),82.2,353.0,8.9e-169,594.3,sp|P82471|GNAQ_RAT Guanine nucleotide-binding protein G(q) subunit alpha OS=Rattus norvegicus OX=10116 GN=Gnaq PE=2 SV=2,GO:0001501; GO:0001508; GO:0001664; GO:0003924; GO:0003925; GO:0005096; GO:0005525; GO:0005737; GO:0005794; GO:0005834; GO:0005886; GO:0005901; GO:0007186; GO:0007189; GO:0007200; GO:0007206; GO:0007207; GO:0007208; GO:0007215; GO:0007218; GO:0007507; GO:0008217; GO:0009791; GO:0010543; GO:0016020; GO:0016322; GO:0021884; GO:0030234; GO:0030425; GO:0031683; GO:0031965; GO:0032024; GO:0034695; GO:0042711; GO:0042733; GO:0043066; GO:0043267; GO:0044297; GO:0044877; GO:0045634; GO:0046872; GO:0047391; GO:0048066; GO:0050821; GO:0060158; GO:0060828; GO:0086100; GO:0098793; GO:0099524; GO:1904888; GO:1990806,action potential [GO:0001508]; adenylate cyclase-activating G protein-coupled receptor signaling pathway [GO:0007189]; cranial skeletal system development [GO:1904888]; developmental pigmentation [GO:0048066]; embryonic digit morphogenesis [GO:0042733]; endothelin receptor signaling pathway [GO:0086100]; forebrain neuron development [GO:0021884]; G protein-coupled receptor signaling pathway [GO:0007186]; glutamate receptor signaling pathway [GO:0007215]; heart development [GO:0007507]; ligand-gated ion channel signaling pathway [GO:1990806]; maternal behavior [GO:0042711]; negative regulation of apoptotic process [GO:0043066]; negative regulation of potassium ion transport [GO:0043267]; neuron remodeling [GO:0016322]; neuropeptide signaling pathway [GO:0007218]; phospholipase C-activating dopamine receptor signaling pathway [GO:0060158]; phospholipase C-activating G protein-coupled acetylcholine receptor signaling pathway [GO:0007207]; phospholipase C-activating G protein-coupled glutamate receptor signaling pathway [GO:0007206]; phospholipase C-activating G protein-coupled receptor signaling pathway [GO:0007200]; phospholipase C-activating serotonin receptor signaling pathway [GO:0007208]; positive regulation of insulin secretion [GO:0032024]; post-embryonic development [GO:0009791]; protein stabilization [GO:0050821]; regulation of blood pressure [GO:0008217]; regulation of canonical Wnt signaling pathway [GO:0060828]; regulation of melanocyte differentiation [GO:0045634]; regulation of platelet activation [GO:0010543]; response to prostaglandin E [GO:0034695]; skeletal system development [GO:0001501],caveola [GO:0005901]; cell body [GO:0044297]; cytoplasm [GO:0005737]; dendrite [GO:0030425]; Golgi apparatus [GO:0005794]; heterotrimeric G-protein complex [GO:0005834]; membrane [GO:0016020]; nuclear membrane [GO:0031965]; plasma membrane [GO:0005886]; postsynaptic cytosol [GO:0099524]; presynapse [GO:0098793],alkylglycerophosphoethanolamine phosphodiesterase activity [GO:0047391]; enzyme regulator activity [GO:0030234]; G protein activity [GO:0003925]; G protein-coupled receptor binding [GO:0001664]; G-protein beta/gamma-subunit complex binding [GO:0031683]; GTP binding [GO:0005525]; GTPase activator activity [GO:0005096]; GTPase activity [GO:0003924]; metal ion binding [GO:0046872]; protein-containing complex binding [GO:0044877],GO:0003824; GO:0003924; GO:0005794; GO:0005829; GO:0005886; GO:0016787; GO:0022414; GO:0043226; GO:0048856; GO:0098772,catalytic activity; GTPase activity; Golgi apparatus; cytosol; plasma membrane; hydrolase activity; reproductive process; organelle; anatomical structure development; molecular function regulator activity
    OG_00248,FUN_001011-T1,Peve_00041060,Pocillopora_meandrina_HIv1___RNAseq.g29164.t1,three_way,70.763,FUN_001011-T1,Q6ZUX7,LHPL2_HUMAN,,LHFPL tetraspan subfamily member 2 protein (Lipoma HMGIC fusion partner-like 2 protein),Homo sapiens (Human),43.1,232.0,6.4e-43,175.6,sp|Q6ZUX7|LHPL2_HUMAN LHFPL tetraspan subfamily member 2 protein OS=Homo sapiens OX=9606 GN=LHFPL2 PE=1 SV=2,GO:0005886; GO:0007338; GO:0016020; GO:0031092; GO:0046545; GO:0046546; GO:1905516,development of primary female sexual characteristics [GO:0046545]; development of primary male sexual characteristics [GO:0046546]; positive regulation of fertilization [GO:1905516]; single fertilization [GO:0007338],membrane [GO:0016020]; plasma membrane [GO:0005886]; platelet alpha granule membrane [GO:0031092],,GO:0005886; GO:0022414,plasma membrane; reproductive process

``` bash
# Count orthologs by each GO category
echo "=== Glycolysis (GO:0006096) ==="
grep -c 'GO:0006096' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Gluconeogenesis (GO:0006094) ==="
grep -c 'GO:0006094' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Lipolysis/lipid catabolism (GO:0016042) ==="
grep -c 'GO:0016042' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Fatty acid beta oxidation (GO:0006635) ==="
grep -c 'GO:0006635' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Starvation (GO:0042594) ==="
grep -c 'GO:0042594' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Lipid biosynthesis (GO:0008610) ==="
grep -c 'GO:0008610' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Protein catabolic process (GO:0030163) ==="
grep -c 'GO:0030163' ../output/13-ortho-select/selected_orthologs.csv || echo "0"

echo "=== Reproductive Process (GO:0022414) ==="
grep -c 'GO:0022414' ../output/13-ortho-select/selected_orthologs.csv || echo "0"
```

    === Glycolysis (GO:0006096) ===
    10
    === Gluconeogenesis (GO:0006094) ===
    13
    === Lipolysis/lipid catabolism (GO:0016042) ===
    21
    === Fatty acid beta oxidation (GO:0006635) ===
    22
    === Starvation (GO:0042594) ===
    12
    === Lipid biosynthesis (GO:0008610) ===
    7
    === Protein catabolic process (GO:0030163) ===
    235
    === Reproductive Process (GO:0022414) ===
    496
