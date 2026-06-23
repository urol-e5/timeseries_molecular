#!/usr/bin/env python3
"""
35-physiology-gene-sets.py

INPUTS
  - 12-ortho-annot annotation_full_go.tsv   (ortholog-group representative -> reviewed Swiss-Prot + GO)
  - 11.3-ortholog-annotation ortholog_annotations_database.csv (OG group -> apul, peve, ptua members)

OUTPUTS (M-multi-species/output/35-physiology-gene-sets/, or argv[1])
  - respiration_geneset_expanded.csv , antioxidant_geneset_expanded.csv ,
    symbiosis_geneset_expanded.csv , growth_geneset_expanded.csv
  Columns: set, group_id, apul, peve, ptua, swissprot_accession, protein_name, evalue, source,
  go_terms, marker_family, subcategory, confidence, justification, reference.
  (biomineralization is built separately by 31.1-biomin-blast.)
"""
import csv, re, os, sys

BASE  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # M-multi-species
ANNOT = os.path.join(BASE, "output/12-ortho-annot/run_20250831_172744/annotation_full_go.tsv")
ORTHO = os.path.join(BASE, "output/11.3-ortholog-annotation/ortholog_annotations_database.csv")
OUT   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "output/35-physiology-gene-sets")
os.makedirs(OUT, exist_ok=True)
EVMAX = 1e-10

REFS = {
 "GO":        "Gene Ontology Consortium 2023, Genetics 224:iyad031 (term membership from reviewed Swiss-Prot annotation)",
 "KEGG_ox":   "KEGG oxidative phosphorylation map00190; Kanehisa & Goto 2000, Nucleic Acids Res 28:27-30",
 "KEGG_tca":  "KEGG citrate (TCA) cycle map00020; Kanehisa & Goto 2000, Nucleic Acids Res 28:27-30",
 "Lesser2006":"Lesser 2006, Annu Rev Physiol 68:253-278 (antioxidant enzyme network in marine invertebrates)",
 "Downs2002": "Downs et al. 2002, Free Radic Biol Med 33:533-543 (coral antioxidant bleaching biomarkers)",
 "Weis2008":  "Weis 2008, J Exp Biol 211:3059-3066 (ROS scavenging in cnidarian bleaching)",
 "Davy2012":  "Davy, Allemand & Weis 2012, Microbiol Mol Biol Rev 76:229-261 (cell biology of cnidarian-dinoflagellate symbiosis)",
 "Pernice2012":"Pernice et al. 2012, ISME J 6:1314-1324 (ammonium assimilation in the symbiosis)",
 "Cui2019":   "Cui et al. 2019, PLoS Genet 15:e1008189 (host nitrogen assimilation/recycling controls symbionts)",
 "Lehnert2014":"Lehnert et al. 2014, G3 4:277-295 (symbiotic vs aposymbiotic transcriptomes)",
 "Hambleton2019":"Hambleton et al. 2019, eLife 8:e43923 (NPC2 sterol transfer in symbiosis)",
 "Dani2017":  "Dani et al. 2017, Cell Microbiol 19:e12753 (NPC1/NPC2 sterol transporters in symbiosis)",
 "Barott2015":"Barott et al. 2015, PNAS 112:607-612 (V-ATPase symbiosome acidification / CCM)",
 "Zoccola2015":"Zoccola et al. 2015, Sci Rep 5:9983 (bicarbonate transporters supply DIC in coral symbiosis)",
 "Bertucci2013":"Bertucci et al. 2013, Bioorg Med Chem 21:1437-1450 (carbonic anhydrases; DIC supply / CCM)",
 "Sproles2018":"Sproles et al. 2018, Mol Phylogenet Evol 120:307-320 (transporter proteins in the symbiosis)",
 "WoodCharlson2006":"Wood-Charlson et al. 2006, Cell Microbiol 8:1985-1993 (lectin/glycan recognition)",
 "Kvennefors2008":"Kvennefors et al. 2008, Dev Comp Immunol 32:1582-1592 (mannose-binding lectin in Acropora)",
 "Chen2003_05":"Chen et al. 2003 BBRC 308:586; 2004 BBRC 324:1024; 2005 BBRC 338:1607 (symbiosome Rab GTPases)",
 "Warner1999": "Warner 1999, Trends Biochem Sci 24:437-440 (economics of ribosome biosynthesis; ribosome content scales with growth/biomass)",
 "Elser2003":  "Elser et al. 2003, Ecol Lett 6:936-943 (growth rate hypothesis: ribosome/rRNA content couples to biomass growth rate)",
 "Saxton2017": "Saxton & Sabatini 2017, Cell 168:960-976 (mTOR signalling in growth and metabolism)",
 "Sonenberg2009":"Sonenberg & Hinnebusch 2009, Cell 136:731-745 (regulation of eukaryotic translation initiation)",
 "Malumbres2009":"Malumbres & Barbacid 2009, Nat Rev Cancer 9:153-166 (cell cycle, cyclins and CDKs)",
}

# ---- (A) domain-defining GO terms: term -> (label, confidence) ----
GO_RESP = {
 "GO:0045333":("cellular respiration","high"), "GO:0009060":("aerobic respiration","high"),
 "GO:0006119":("oxidative phosphorylation","high"), "GO:0006099":("tricarboxylic acid cycle","high"),
 "GO:0022904":("respiratory electron transport chain","high"),
 "GO:0006120":("mito ET, NADH to ubiquinone","high"), "GO:0006121":("mito ET, succinate to ubiquinone","high"),
 "GO:0006122":("mito ET, ubiquinol to cyt c","high"), "GO:0006123":("mito ET, cyt c to oxygen","high"),
 "GO:0042775":("mito ATP synth coupled ET","high"), "GO:0008137":("NADH dehydrogenase (ubiquinone) activity","high"),
 "GO:0004129":("cytochrome-c oxidase activity","high"), "GO:0046933":("proton-transporting ATP synthase activity","high"),
 "GO:0015980":("energy derivation by oxidation","moderate"),
}
GO_ANTI = {
 "GO:0016209":("antioxidant activity","high"), "GO:0004601":("peroxidase activity","high"),
 "GO:0004784":("superoxide dismutase activity","high"), "GO:0042744":("hydrogen peroxide catabolic process","high"),
 "GO:0006749":("glutathione metabolic process","high"), "GO:0045454":("cell redox homeostasis","high"),
 "GO:0004362":("glutathione-disulfide reductase activity","high"), "GO:0004602":("glutathione peroxidase activity","high"),
 "GO:0098869":("cellular oxidant detoxification","high"), "GO:0008379":("thioredoxin peroxidase activity","high"),
 "GO:0032928":("regulation of superoxide metabolic process","moderate"),
}
GO_SYMB = {  # sparse; supporting only
 "GO:0072488":("ammonium transmembrane transport","high"), "GO:0015696":("ammonium transport","high"),
}
GO_GROWTH = {  # anabolic / protein-biosynthetic capacity = molecular correlate of biomass + protein
 "GO:0042254":("ribosome biogenesis","high"), "GO:0006364":("rRNA processing","high"),
 "GO:0003735":("structural constituent of ribosome","high"), "GO:0002181":("cytoplasmic translation","high"),
 "GO:0006412":("translation","high"), "GO:0006413":("translational initiation","high"),
 "GO:0006414":("translational elongation","high"), "GO:0006418":("tRNA aminoacylation for translation","high"),
 "GO:0004812":("aminoacyl-tRNA ligase activity","high"), "GO:0031929":("TOR signaling","high"),
}

# ---- (B) curated families: (label, regex, subcat, confidence, justification, refkey) ----
CUR_RESP = [
 ("Complex I (NADH:ubiquinone oxidoreductase)", r"NADH dehydrogenase \[ubiquinone\]|NADH-ubiquinone oxidoreductase",
  "ETC complex I","high","Core subunit of respiratory complex I.","KEGG_ox"),
 ("Complex II (succinate dehydrogenase)", r"Succinate dehydrogenase \[ubiquinone\]",
  "ETC complex II","high","Succinate dehydrogenase couples TCA to ETC.","KEGG_ox"),
 ("Complex III (cytochrome bc1)", r"Cytochrome b-c1 complex subunit|Cytochrome c1, heme",
  "ETC complex III","high","Cytochrome bc1 subunit of the respiratory chain.","KEGG_ox"),
 ("Complex IV (cytochrome c oxidase)", r"Cytochrome c oxidase subunit|Cytochrome c oxidase polypeptide",
  "ETC complex IV","high","Terminal oxidase reducing O2 to water.","KEGG_ox"),
 ("Cytochrome c", r"^Cytochrome c$|^Cytochrome c,|^Cytochrome c-",
  "ETC electron carrier","high","Soluble carrier between complexes III and IV.","KEGG_ox"),
 ("ATP synthase (complex V)", r"ATP synthase subunit|ATP synthase F\(?[01]\)?|ATP synthase peripheral stalk|ATP synthase-coupling",
  "ETC complex V","high","F1Fo-ATP synthase subunit; synthesises ATP from the proton-motive force.","KEGG_ox"),
 ("TCA cycle enzymes", r"^Citrate synthase(?!-lysine)|aconitate hydratase, mitochondrial|Isocitrate dehydrogenase|2-oxoglutarate dehydrogenase|oxoglutarate dehydrogenase complex|Succinate--CoA ligase|Fumarate hydratase|Malate dehydrogenase, mitochondrial|Dihydrolipoyllysine-residue succinyltransferase component of 2-oxoglutarate",
  "TCA cycle","high","Core TCA-cycle enzyme generating reducing equivalents for the ETC.","KEGG_tca"),
 ("Pyruvate dehydrogenase complex", r"Pyruvate dehydrogenase|Dihydrolipoyllysine-residue acetyltransferase|Dihydrolipoyl dehydrogenase",
  "TCA feeder","high","PDH complex converts pyruvate to acetyl-CoA, feeding the TCA cycle / aerobic respiration.","KEGG_tca"),
 ("ETF-ubiquinone oxidoreductase", r"Electron transfer flavoprotein-ubiquinone oxidoreductase",
  "ETC electron carrier","moderate","Feeds electrons from flavoprotein dehydrogenases to ubiquinone.","KEGG_ox"),
 ("ADP/ATP translocase", r"ADP/ATP translocase|ADP,ATP carrier",
  "OXPHOS support","moderate","Exchanges mitochondrial ATP for cytosolic ADP, sustaining OXPHOS flux.","KEGG_ox"),
]
CUR_ANTI = [
 ("Superoxide dismutase", r"Superoxide dismutase \[", "ROS detoxification","high",
  "SOD dismutes superoxide to H2O2; front-line antioxidant.","Downs2002"),
 ("Catalase", r"^Catalase(?! ?-)", "ROS detoxification","high",
  "Catalase decomposes H2O2; core antioxidant.","Downs2002"),
 ("Glutathione peroxidase", r"Glutathione peroxidase", "ROS detoxification","high",
  "GPx reduces peroxides using glutathione.","Lesser2006"),
 ("Peroxiredoxin", r"Peroxiredoxin|Thioredoxin-dependent peroxide reductase", "ROS detoxification","high",
  "Thioredoxin-linked peroxidase, abundant under oxidative stress.","Lesser2006"),
 ("Thioredoxin / reductase", r"^Thioredoxin(?!-related transmembrane| domain-containing|-like)|Thioredoxin reductase",
  "thiol-redox","moderate","Maintains thiol-redox balance; recycles peroxiredoxins.","Lesser2006"),
 ("Glutaredoxin", r"Glutaredoxin", "thiol-redox","moderate","Glutathione-dependent thiol oxidoreductase.","Lesser2006"),
 ("Glutathione reductase", r"Glutathione reductase", "glutathione system","high",
  "Regenerates reduced glutathione, sustaining antioxidant capacity.","Lesser2006"),
 ("Glutathione synthesis", r"Glutathione synthetase|Glutamate--cysteine ligase", "glutathione system","moderate",
  "Sets the cellular GSH pool underlying total antioxidant capacity.","Lesser2006"),
 ("Glutathione S-transferase", r"[Gg]lutathione S-transferase", "conjugation/detox","moderate",
  "Conjugates GSH to oxidised products.","Lesser2006"),
 ("Peroxidase (other)", r"Peroxidasin|Glutathione peroxidase|Thioredoxin peroxidase|Ascorbate peroxidase",
  "ROS detoxification","moderate","Additional peroxidase scavenging activity.","Lesser2006"),
 ("Heme oxygenase", r"Heme oxygenase", "oxidative-stress response","moderate",
  "HO-1 is a cytoprotective, oxidative-stress-induced antioxidant enzyme.","Lesser2006"),
 ("Glucose-6-phosphate dehydrogenase", r"Glucose-6-phosphate 1-dehydrogenase", "NADPH supply","moderate",
  "Supplies NADPH that regenerates glutathione/thioredoxin antioxidant systems.","Lesser2006"),
 ("Methionine sulfoxide reductase", r"Methionine sulfoxide reductase|Peptide methionine sulfoxide reductase",
  "oxidative-damage repair","moderate","Repairs oxidised methionine residues.","Lesser2006"),
 ("Sulfiredoxin", r"Sulfiredoxin", "thiol-redox","moderate","Regenerates hyperoxidised peroxiredoxins.","Lesser2006"),
 ("Ferritin", r"Ferritin", "iron sequestration","moderate",
  "Sequesters iron, limiting Fenton-reaction ROS generation.","Lesser2006"),
 ("Superoxide dismutase copper chaperone", r"Copper chaperone for superoxide dismutase",
  "ROS detoxification","moderate","Delivers copper to Cu/Zn-SOD.","Lesser2006"),
]
CUR_SYMB = [
 ("Ammonium transporter / Rhesus", r"Ammonium transporter|Rhesus|Rh type", "N exchange","high",
  "Host ammonium/Rh transporters mediate symbiotic nitrogen exchange.","Pernice2012"),
 ("Glutamine synthetase", r"Glutamine synthetase|Glutamate--ammonia ligase", "N assimilation","high",
  "GS assimilates ammonium into glutamine; central to host control of symbiont N.","Cui2019"),
 ("Glutamate dehydrogenase", r"Glutamate dehydrogenase", "N assimilation","high",
  "GDH interconverts glutamate/2-oxoglutarate+ammonium in symbiotic N recycling.","Cui2019"),
 ("Glutamate synthase (GOGAT)", r"Glutamate synthase", "N assimilation","high",
  "GOGAT completes the GS/GOGAT ammonium-assimilation cycle.","Cui2019"),
 ("NPC1/NPC2 sterol transporter", r"NPC intracellular cholesterol transporter|Niemann-Pick", "sterol transfer","high",
  "Niemann-Pick C proteins mediate host-symbiont sterol transfer.","Hambleton2019"),
 ("Oxysterol-binding protein", r"Oxysterol-binding", "sterol transfer","moderate",
  "Sterol-binding/transfer proteins implicated in symbiotic lipid trafficking.","Dani2017"),
 ("Carbonic anhydrase", r"[Cc]arbonic anhydrase", "DIC supply / CCM","moderate",
  "Interconverts CO2/HCO3- to supply inorganic carbon for symbiont photosynthesis (overlaps biomin).","Bertucci2013"),
 ("Bicarbonate / anion transporter", r"Anion exchange protein|Sodium bicarbonate|Sodium-driven chloride bicarbonate|Electrogenic sodium bicarbonate|Solute carrier family 4 member",
  "DIC supply / CCM","moderate","Bicarbonate transporters (SLC4 family) deliver DIC to the symbiont.","Zoccola2015"),
 ("V-type proton ATPase", r"V-type proton ATPase", "symbiosome CCM","moderate",
  "Acidifies the symbiosome to concentrate carbon for photosynthesis.","Barott2015"),
 ("Facilitated glucose transporter", r"facilitated glucose transporter", "photosynthate transfer","moderate",
  "Candidate conduit for translocated photosynthate (glucose).","Sproles2018"),
 ("Aquaporin / aquaglyceroporin", r"Aquaporin", "nutrient/water transport","moderate",
  "Aquaglyceroporins implicated in glycerol/water/DIC exchange in the symbiosis.","Sproles2018"),
 ("Mannose / C-type lectin receptor", r"Macrophage mannose receptor|C-type mannose receptor|C-type lectin domain family",
  "symbiont recognition","moderate","Host lectins recognise symbiont surface glycans.","WoodCharlson2006"),
 ("Tachylectin / ficolin / fucolectin", r"Tachylectin|Ficolin|Fucolectin", "symbiont recognition","moderate",
  "Carbohydrate-recognition lectins implicated in symbiont uptake/winnowing.","Kvennefors2008"),
 ("Symbiosome Rab GTPase (Rab5/7/11)", r"Ras-related protein Rab-5\b|Ras-related protein Rab-7\b|Ras-related protein Rab-11",
  "symbiosome trafficking","moderate","Rab5/7/11 localise to the symbiosome and regulate its biogenesis.","Chen2003_05"),
]
# growth/biomass = anabolic capacity. TOR family is listed FIRST so 'Ribosomal protein S6
# kinase' matches growth-signalling, not the generic ribosomal-protein family.
CUR_GROWTH = [
 ("TOR pathway (growth signalling)", r"Serine/threonine-protein kinase (mTOR|TOR)|Target of rapamycin|Regulatory-associated protein of mTOR|Rapamycin-insensitive companion|Ribosomal protein S6 kinase|Eukaryotic translation initiation factor 4E-binding",
  "growth signalling (TOR)","high","TOR-pathway regulator linking nutrient status to ribosome biogenesis and cap-dependent translation; master controller of cell growth/anabolism.","Saxton2017"),
 ("Ribosomal protein (cytoplasmic)", r"[Rr]ibosomal (?:subunit )?protein(?! S6 kinase)",
  "translation (ribosome)","high","Cytoplasmic ribosomal protein; ribosome content is the canonical molecular correlate of biomass/protein growth rate (growth-rate hypothesis).","Warner1999"),
 ("Ribosome biogenesis / rRNA processing", r"[Rr]ibosome biogenesis|Pre-rRNA-processing|rRNA[ -]processing|Ribosomal RNA processing|Ribosome assembly|Ribosome maturation|Ribosome production|Ribosome-binding factor",
  "ribosome biogenesis","high","Ribosome-biogenesis / rRNA-processing factor; ribosome production sets anabolic capacity.","Warner1999"),
 ("Translation initiation factor", r"Eukaryotic translation initiation factor|Eukaryotic initiation factor",
  "translation (initiation)","high","Cytoplasmic translation initiation factor (eIF); rate-controlling for protein synthesis.","Sonenberg2009"),
 ("Translation elongation factor", r"Elongation factor 1-|Elongation factor 2|Eukaryotic translation elongation|Eukaryotic peptide chain release",
  "translation (elongation)","high","Cytoplasmic translation elongation / release factor (eEF1/eEF2).","Sonenberg2009"),
 ("Aminoacyl-tRNA synthetase", r"--tRNA ligase|tRNA synthetase|Aminoacyl tRNA synthase",
  "translation (tRNA charging)","high","Aminoacyl-tRNA synthetase charges tRNAs for protein synthesis.","GO"),
 ("Replication licensing (MCM) / PCNA", r"Proliferating cell nuclear antigen|DNA replication licensing factor MCM|Minichromosome maintenance",
  "cell proliferation","moderate","DNA-replication / proliferation marker (PCNA, MCM2-7); proxy for the cell-division contribution to tissue biomass.","Malumbres2009"),
 ("Cell-cycle cyclin", r"^Cyclin-[ABDE]\b|G1/S-specific cyclin|G2/mitotic-specific cyclin",
  "cell proliferation","moderate","Cell-cycle (A/B/D/E) cyclin driving cell-division progression.","Malumbres2009"),
 ("Cell-cycle CDK", r"Cyclin-dependent kinase [1246]\b",
  "cell proliferation","moderate","Cell-cycle cyclin-dependent kinase (CDK1/2/4/6).","Malumbres2009"),
]

# ---- name-collision vetoes (drop matches that are NOT the intended marker) ----
EXCLUDE = {
 "respiration":[(r"assembly factor|assembly protein|chaperone|copper chaperone|intermediate-associated","assembly/chaperone, not catalytic"),
   (r"methyltransferase|hydroxylase NDUFAF|Arginine-hydroxylase|lysine N-methyltransferase","modifier of an ETC subunit"),
   (r"Cytochrome b-245|NADPH oxidase|Superoxide-generating","NADPH oxidase (ROS-generating), not respiration"),
   (r"Cytoplasmic aconitate|Iron regulatory protein|Iron-responsive","cytosolic aconitase/IRP1, not TCA")],
 "antioxidant":[(r"Thioredoxin-related transmembrane|Thioredoxin domain-containing|Thioredoxin-like protein 4|U5 snRNP|disulfide-isomerase|ERp44|spliceosomal","thioredoxin-fold non-antioxidant"),
   (r"Cytoglobin","globin with non-antioxidant roles"),
   (r"Peroxidasin.*(?:motif|domain only)","ECM peroxidasin, ambiguous")],
 "symbiosis":[(r"Rab-7L1|Rab-29","not the symbiosome Rab7"),
   (r"V-type proton ATPase.*(assembly|accessory|S1 accessory)","V-ATPase accessory factor")],
 "growth":[(r"mitochondrial|chloroplastic|Ribosome-recycling factor",
    "mitochondrial/organellar translation apparatus (overlaps respiration; growth set is cytoplasmic biomass)"),
   (r"Transcription (elongation|initiation|termination) factor|General transcription factor|Transcription factor|RNA polymerase",
    "transcription apparatus, not translation"),
   (r"Eukaryotic translation initiation factor 2-alpha kinase|eIF-2-alpha kinase|Serine/threonine-protein kinase GCN2|PERK\b",
    "stress-induced translational repressor (GCN2/PERK/PKR/HRI), not anabolic"),
   (r"Cyclin-dependent kinase inhibitor|Cyclin-dependent kinase [789]\b|Cyclin-dependent kinase 1[0-3]\b|Cyclin-dependent kinase-like|Cyclin-[CHKLTGY]\b|Cyclin-dependent kinases regulatory",
    "CDK inhibitor or transcriptional/non-cell-cycle CDK/cyclin"),
   (r"chaperone|Heat shock|T-complex protein",
    "protein-folding chaperone (stress-confounded), not anabolic biosynthesis")],
}

GO_TERMS = {"respiration":GO_RESP, "antioxidant":GO_ANTI, "symbiosis":GO_SYMB, "growth":GO_GROWTH}
CURATED  = {"respiration":CUR_RESP, "antioxidant":CUR_ANTI, "symbiosis":CUR_SYMB, "growth":CUR_GROWTH}


def parse_ev(s):
    try: return float(s)
    except: return 1.0


def present(x):
    """True if an ortholog-member id is real (non-empty, not NA)."""
    return bool(x and x.strip() and x.strip().upper() != "NA")


def load_ortho():
    """Apul FUN id -> (group_id, peve, ptua)."""
    m = {}
    for row in csv.DictReader(open(ORTHO)):
        a = (row.get("apul") or "").strip()
        if a and a.upper() != "NA":
            m[a] = (row.get("group_id",""), row.get("peve",""), row.get("ptua",""))
    return m


def build():
    omap = load_ortho()
    annot = list(csv.DictReader(open(ANNOT), delimiter="\t"))
    summary = {}
    for dm in ["respiration", "antioxidant", "symbiosis", "growth"]:
        gos, cur, exc = GO_TERMS[dm], CURATED[dm], EXCLUDE[dm]
        out = {}
        for r in annot:
            name = (r.get("protein_name") or "").strip()
            if not name:
                continue
            if any(re.search(p, name) for p, _ in exc):        # collision veto
                continue
            if r.get("Reviewed", "").strip().lower() != "reviewed":
                continue
            if parse_ev(r.get("evalue", "1")) > EVMAX:
                continue
            gids = set((r.get("go_ids", "") or "").replace(" ", "").split(";"))
            go_hits = [(t, gos[t][0], gos[t][1]) for t in gos if t in gids]
            fam = None
            for lab, rx, sub, conf, just, ref in cur:
                if re.search(rx, name):
                    fam = (lab, sub, conf, just, ref); break
            if not go_hits and not fam:                         # needs source A or B
                continue
            fun = (r.get("query") or "").strip()
            grp, peve, ptua = omap.get(fun, ("", "", ""))
            # OPTION B (species-symmetric): keep ONLY three-way orthologs -- groups present in
            # all three species -- so every species' set is the IDENTICAL group set and the
            # cross-species comparison treats Apul/Peve/Ptua on equal footing. Output is keyed on
            # the ortholog group, with symmetric apul/peve/ptua member columns (no Apul "baseline").
            # NB: the functional annotation (protein_name/GO) is the ortholog-group representative's
            # -- the Apul member for three-way groups; orthology conserves function, so the call
            # applies to the whole group. (Fully independent per-species annotation = Option A.)
            if not (present(grp) and present(peve) and present(ptua)):
                continue
            if grp in out:                                      # one row per ortholog group
                continue
            confs = [c for _, _, c in go_hits] + ([fam[2]] if fam else [])
            conf = "high" if "high" in confs else "moderate"
            source = "GO+curated" if (go_hits and fam) else ("GO" if go_hits else "curated")
            go_str = "; ".join(sorted({f"{t}({lab})" for t, lab, _ in go_hits}))
            if fam:
                sub, just, ref, family = fam[1], fam[3], fam[4], fam[0]
            else:                                               # GO-only group
                sub = go_hits[0][1]; family = f"GO:{go_hits[0][1]}"
                just = f"Annotated with {go_str}."; ref = "GO"
            out[grp] = dict(set=dm, group_id=grp, apul=fun, peve=peve, ptua=ptua,
                swissprot_accession=r.get("accession", ""), protein_name=name, evalue=r.get("evalue", ""),
                source=source, go_terms=go_str, marker_family=family, subcategory=sub,
                confidence=conf, justification=just, reference=REFS.get(ref, ref))
        rows = sorted(out.values(),
                      key=lambda d: (d["confidence"] != "high", d["subcategory"], d["marker_family"], d["group_id"]))
        cols = ["set","group_id","apul","peve","ptua","swissprot_accession","protein_name","evalue",
                "source","go_terms","marker_family","subcategory","confidence","justification","reference"]
        with open(os.path.join(OUT, f"{dm}_geneset_expanded.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
        nh = sum(r["confidence"] == "high" for r in rows)
        bysrc = {s: sum(r["source"] == s for r in rows) for s in ("GO", "curated", "GO+curated")}
        summary[dm] = (len(rows), nh, bysrc)
        print(f"{dm:12s}: {len(rows):3d} three-way ortholog groups  high={nh}  "
              f"source GO={bysrc['GO']} curated={bysrc['curated']} both={bysrc['GO+curated']}")
    return summary


if __name__ == "__main__":
    build()
