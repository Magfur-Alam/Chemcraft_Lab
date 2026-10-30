import io, json, textwrap, requests, streamlit as st, pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from streamlit.components.v1 import html

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
st.set_page_config(page_title="ChemInformatics Virtual Laboratory", page_icon="🧪", layout="wide")

# ---------- Theme ----------
def apply_theme():
    themes = {
        "Midnight Blue": ("#0b1220", "#16243a", "#53b7ff", "#f4f7fb"),
        "Emerald": ("#071b18", "#12352d", "#36d6a0", "#f2fff9"),
        "Purple": ("#171126", "#2b1d45", "#bd91ff", "#fbf7ff"),
        "Light": ("#f4f7fb", "#ffffff", "#1769aa", "#152235"),
        "Rose": ("#24131b", "#3a202d", "#ff8db2", "#fff5f8"),
    }
    bg, panel, accent, fg = themes[st.session_state.get("theme", "Midnight Blue")]
    st.markdown(f"""
    <style>
    .stApp {{background:{bg}; color:{fg}}}
    [data-testid="stSidebar"] {{background:{panel}}}
    h1,h2,h3 {{color:{accent}}}
    .stButton>button {{border-radius:10px; border:1px solid {accent};}}
    div[data-testid="stMetric"] {{background:{panel}; padding:12px; border-radius:12px;}}
    </style>""", unsafe_allow_html=True)

if "theme" not in st.session_state: st.session_state.theme = "Midnight Blue"
apply_theme()

# ---------- Chemistry helpers ----------
@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_search(term, mode="name"):
    term = term.strip()
    if not term: return None
    safe = requests.utils.quote(term, safe="")
    url = f"{PUBCHEM}/compound/{mode}/{safe}/cids/JSON"
    r = requests.get(url, timeout=25)
    if r.status_code == 404: return []
    r.raise_for_status()
    return r.json().get("IdentifierList", {}).get("CID", [])

@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_properties(cid):
    props = "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,IUPACName,InChIKey"
    url = f"{PUBCHEM}/compound/cid/{int(cid)}/property/{props}/JSON"
    r = requests.get(url, timeout=25); r.raise_for_status()
    return r.json()["PropertyTable"]["Properties"][0]

@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_sdf(cid, record_type="3d"):
    url = f"{PUBCHEM}/compound/cid/{int(cid)}/record/SDF"
    r = requests.get(url, params={"record_type": record_type}, timeout=30)
    if r.status_code != 200: return ""
    return r.text

def molecule_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES. Please check the notation.")
    return mol

def descriptors(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "Exact mass": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond donors": rdMolDescriptors.CalcNumHBD(mol),
        "H-bond acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "Rotatable bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "Heavy atoms": mol.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(mol),
    }

def image_2d(mol, size=(650, 430)):
    m = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(m)
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    drawer.DrawMolecule(m); drawer.FinishDrawing()
    return drawer.GetDrawingText()

def make_3d(mol):
    m = Chem.AddHs(Chem.Mol(mol))
    status = AllChem.EmbedMolecule(m, randomSeed=42)
    if status != 0: return None
    AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    return Chem.MolToMolBlock(m)

def show_3d(molblock, style="stick", colorscheme="Jmol", background="#101522"):
    # 3Dmol.js is loaded in the browser; the molecule is generated locally by RDKit.
    sdf = json.dumps(molblock)
    html_code = f"""
    <div id="viewer" style="width:100%;height:470px;position:relative;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    const element=document.getElementById('viewer');
    const viewer=$3Dmol.createViewer(element,{{backgroundColor:'{background}'}});
    const data={sdf};
    viewer.addModel(data,'mol');
    viewer.setStyle({{}},{{{style}:{{colorscheme:'{colorscheme}'}}}});
    viewer.zoomTo(); viewer.render(); viewer.zoom(0.9);
    </script>"""
    html(html_code, height=490, scrolling=False)

def report_pdf(title, smiles, props, cid=None):
    from reportlab.lib.pagesizes import A4
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib import colors
    buf=io.BytesIO(); doc=SimpleDocTemplate(buf, pagesize=A4)
    styles=getSampleStyleSheet(); story=[Paragraph("ChemInformatics Virtual Laboratory",styles["Title"]),
        Paragraph(f"Molecular Report: {title}",styles["Heading2"]), Spacer(1,12),
        Paragraph(f"<b>SMILES:</b> {smiles}",styles["BodyText"]), Spacer(1,12)]
    if cid: story.append(Paragraph(f"<b>PubChem CID:</b> {cid}",styles["BodyText"]))
    data=[["Property","Value"]]+[[str(k),str(v)] for k,v in props.items()]
    table=Table(data, colWidths=[230,260]); table.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#163b5c")),("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("GRID",(0,0),(-1,-1),0.5,colors.grey),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#eef4f8")]),
        ("PADDING",(0,0),(-1,-1),7)]))
    story.append(table); story.append(Spacer(1,16))
    story.append(Paragraph("Calculated locally with RDKit. PubChem data, when shown, is retrieved from the PubChem PUG REST service.",styles["Italic"]))
    doc.build(story); buf.seek(0); return buf.getvalue()

# ---------- Sidebar ----------
with st.sidebar:
    st.markdown("# 🧪 Virtual Lab")
    page = st.radio("Navigation", ["🏠 Home","📖 Theory","🔎 PubChem Search","🧬 Molecular Visualization",
        "📊 Molecular Properties Report","📈 Structure–Property Analysis","📝 Assessment","📄 Final Project Report"])
    st.divider()
    st.subheader("🎓 Student Information")
    student = st.text_input("Student Name", key="student")
    regno = st.text_input("Registration Number", key="regno")
    st.divider()
    st.subheader("🎨 Appearance")
    st.selectbox("Choose colour theme", ["Midnight Blue","Emerald","Purple","Light","Rose"],
                 key="theme", on_change=apply_theme)
    st.divider()
    st.subheader("🎓 Learning Objectives")
    st.write("Students can learn:")
    for objective in ["SMILES notation", "Molecular representation", "2D visualization",
                      "3D visualization", "Molecular descriptors",
                      "Structure–property relationships"]:
        st.write("• " + objective)
    st.caption("Data source: PubChem • Chemistry calculations: RDKit")

# Store selected molecule between pages
if "smiles" not in st.session_state: st.session_state.smiles = "CCO"
if "mol_name" not in st.session_state: st.session_state.mol_name = "Ethanol"
if "cid" not in st.session_state: st.session_state.cid = None

def set_molecule(name, smiles, cid=None):
    st.session_state.mol_name=name; st.session_state.smiles=smiles; st.session_state.cid=cid

def current_mol():
    return molecule_from_smiles(st.session_state.smiles)

# ---------- Pages ----------
if page == "🏠 Home":
    st.title("🧪 ChemInformatics Virtual Laboratory")
    st.subheader("SRM Institute of Science and Technology, Tiruchirappalli")
    st.write("An interactive learning workspace for SMILES, molecular structures, chemical descriptors and PubChem compound records.")
    c1,c2,c3=st.columns(3)
    c1.metric("Structure tools","2D + 3D")
    c2.metric("Database","PubChem PUG REST")
    c3.metric("Chemistry engine","RDKit")
    st.markdown("### What can students do?")
    for item in ["Search PubChem by compound name or CID","Enter and validate SMILES notation",
        "Generate 2D structures and interactive 3D models","Explore molecular properties and descriptors",
        "Compare molecules and plot property relationships","Download molecular reports and complete quizzes"]:
        st.write("✅ "+item)
    st.info("Try the example compounds below or open PubChem Search from the sidebar.")
    examples={"Ethanol":"CCO","Benzene":"c1ccccc1","Phenol":"Oc1ccccc1","Acetic acid":"CC(=O)O","Aspirin":"CC(=O)Oc1ccccc1C(=O)O","Caffeine":"Cn1c(=O)c2c(ncn2C)n(C)c1=O"}
    st.markdown("### Example Molecules")
    example_rows = [{"Molecule": name, "SMILES": smi} for name, smi in examples.items()]
    st.dataframe(pd.DataFrame(example_rows), hide_index=True, use_container_width=True)
    st.caption("Select a molecule to load it into the visualization and properties tools.")
    cols=st.columns(3)
    for i,(name,smi) in enumerate(examples.items()):
        with cols[i%3]:
            if st.button("Load " + name, key="ex_"+name):
                set_molecule(name,smi); st.rerun()
    st.divider()
    st.markdown("### 👨‍🏫 Developed by")
    st.markdown("**Dr. C. Manivannan**")
    st.write("Assistant Professor")
    st.write("Division of Chemistry, School of Sciences")
    st.write("Faculty of Engineering and Technology")
    st.write("SRM Institute of Science and Technology, Tiruchirappalli")
    st.caption("For academic queries, contact the faculty through official SRM channels.")

elif page == "📖 Theory":
    st.title("📖 Theory & Learning")
    topics={
    "What is Chemoinformatics?":"Chemoinformatics uses computational and information methods to store, search, analyse and visualize chemical information.",
    "SMILES notation":"SMILES is a text representation of a molecular structure. Example: CCO represents ethanol. Atoms are written as symbols; branches use parentheses and rings use matching digits.",
    "2D and 3D structures":"A 2D drawing shows atoms and bonds on a plane. A 3D model adds spatial coordinates and helps explore molecular shape and conformations.",
    "Molecular descriptors":"Descriptors are numerical or categorical values that summarize features such as molecular weight, lipophilicity (LogP), polar surface area and hydrogen-bond counts.",
    "Structure–property relationship":"A structure–property study compares molecular features with properties. A correlation in a small teaching dataset does not prove causation.",
    "PubChem":"PubChem is an open chemistry database from the US National Library of Medicine. Its PUG REST API provides programmatic access to compound identifiers, properties and records."
    }
    for title,body in topics.items():
        with st.expander(title): st.write(body)
    st.caption("Educational content only. Calculated descriptors are not a substitute for experimental measurements.")

elif page == "🔎 PubChem Search":
    st.title("🔎 Search PubChem")
    st.write("Search the PubChem Compound database by common name or PubChem CID.")
    mode=st.radio("Search by",["Compound name","PubChem CID"],horizontal=True)
    term=st.text_input("Enter compound name or CID",placeholder="e.g. aspirin or 2244")
    if st.button("Search PubChem", type="primary"):
        try:
            with st.spinner("Contacting PubChem..."):
                cids=[int(term)] if mode=="PubChem CID" and term.strip().isdigit() else (
                    pubchem_search(term,"name") if mode=="Compound name" else [])
                if not cids: st.warning("No matching CID found. Try another spelling or search by CID.")
                else:
                    st.session_state.search_cids=cids[:20]
                    st.success(f"Found {len(cids)} matching record(s). Showing up to 20.")
        except requests.RequestException as e: st.error(f"PubChem request failed: {e}")
    for cid in st.session_state.get("search_cids",[]):
        try:
            p=pubchem_properties(cid)
            with st.container(border=True):
                st.markdown(f"**{p.get('IUPACName','Compound')}**  •  CID {cid}")
                st.write("Formula:",p.get("MolecularFormula","—"),"| Molecular weight:",p.get("MolecularWeight","—"))
                st.code(p.get("CanonicalSMILES",p.get("IsomericSMILES","")),language="text")
                if st.button(f"Load CID {cid} into laboratory",key=f"load{cid}"):
                    set_molecule(p.get("IUPACName",f"CID {cid}"),p.get("CanonicalSMILES",p.get("IsomericSMILES","")),cid)
                    st.rerun()
                st.markdown(f"[Open PubChem record](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})")
        except Exception as e: st.warning(f"Could not load CID {cid}: {e}")
    st.caption("PubChem is a very large and continuously updated database. Search results depend on PubChem indexing and API limits; this app does not download the entire database.")

elif page == "🧬 Molecular Visualization":
    st.title("🧬 Molecular Visualization")
    mode=st.radio("Input source",["Use current molecule","Enter SMILES"],horizontal=True)
    if mode=="Enter SMILES":
        name=st.text_input("Molecule name","Custom molecule")
        smi=st.text_area("SMILES notation",st.session_state.smiles)
        if st.button("Validate and load"):
            try: molecule_from_smiles(smi); set_molecule(name,smi); st.success("Valid SMILES loaded.")
            except ValueError as e: st.error(str(e))
    try:
        mol=current_mol()
        st.markdown(f"### {st.session_state.mol_name}")
        st.code(st.session_state.smiles)
        left,right=st.columns(2)
        with left:
            st.subheader("2D structure")
            st.image(image_2d(mol), use_container_width=True)
        with right:
            st.subheader("Interactive 3D structure")
            style=st.selectbox("Visualization style",["stick","sphere","line","cross"])
            colorscheme=st.selectbox("Atom colour scheme",["Jmol","Rasmol","cyanCarbon","greenCarbon","magentaCarbon","yellowCarbon","orangeCarbon","whiteCarbon"])
            st.caption("Drag to rotate • scroll/pinch to zoom")
            block=make_3d(mol)
            if block: show_3d(block,style,colorscheme)
            else: st.warning("Could not generate a 3D conformer for this structure.")
        st.success("SMILES validation successful.")
    except Exception as e: st.error(f"Unable to visualize molecule: {e}")

elif page == "📊 Molecular Properties Report":
    st.title("📊 Molecular Properties Report")
    try:
        mol=current_mol(); props=descriptors(mol)
        st.subheader(st.session_state.mol_name)
        st.dataframe(pd.DataFrame(list(props.items()),columns=["Property","Value"]),hide_index=True,use_container_width=True)
        st.download_button("⬇️ Download molecular report (PDF)",
            data=report_pdf(st.session_state.mol_name,st.session_state.smiles,props,st.session_state.cid),
            file_name="molecular_report.pdf",mime="application/pdf")
        if st.session_state.cid: st.markdown(f"[View source record on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{st.session_state.cid})")
    except Exception as e: st.error(str(e))

elif page == "📈 Structure–Property Analysis":
    st.title("📈 Structure–Property Analysis")
    st.write("Compare descriptors for a few example molecules. Select molecules to create a simple teaching plot.")
    examples={"Ethanol":"CCO","Benzene":"c1ccccc1","Phenol":"Oc1ccccc1","Acetic acid":"CC(=O)O","Aspirin":"CC(=O)Oc1ccccc1C(=O)O","Caffeine":"Cn1c(=O)c2c(ncn2C)n(C)c1=O"}
    chosen=st.multiselect("Choose molecules",list(examples),default=["Ethanol","Benzene","Phenol"])
    metric=st.selectbox("Property to plot",["Molecular weight (g/mol)","LogP","TPSA (Å²)","H-bond donors","H-bond acceptors","Rotatable bonds"])
    rows=[]
    for name in chosen:
        d=descriptors(molecule_from_smiles(examples[name])); rows.append({"Molecule":name,metric:d[metric]})
    if rows:
        df=pd.DataFrame(rows); st.dataframe(df,hide_index=True,use_container_width=True)
        fig,ax=plt.subplots(); ax.bar(df["Molecule"],df[metric]); ax.set_ylabel(metric); ax.tick_params(axis="x",rotation=25); fig.tight_layout()
        st.pyplot(fig); st.download_button("Download analysis CSV",df.to_csv(index=False), "structure_property_analysis.csv","text/csv")
    else: st.info("Select at least one molecule.")

elif page == "📝 Assessment":
    st.title("📝 Chemistry & Chemoinformatics Assessment")
    questions=[
      ("What does SMILES represent?",["A text representation of molecular structure","A protein sequence","A unit of pressure"],0),
      ("Which tool is used here for cheminformatics calculations?",["NumPy","RDKit","PowerPoint"],1),
      ("What does molecular weight measure?",["Total mass based on atomic weights","Number of rings only","Molecule colour"],0),
      ("What is PubChem CID?",["A unique compound identifier","A type of chemical bond","A graph axis"],0),
      ("What does TPSA describe?",["Topological polar surface area","Temperature","Total protein sequence"],0)]
    answers=[]
    for i,(q,opts,correct) in enumerate(questions):
        ans=st.radio(f"{i+1}. {q}",opts,key=f"quiz{i}")
        answers.append(opts.index(ans))
    if st.button("Submit assessment"):
        score=sum(a==q[2] for a,q in zip(answers,questions))
        st.success(f"Your score: {score}/{len(questions)}")
        if score==len(questions): st.balloons()
        st.caption("Review the Theory section for topics you want to practise.")

elif page == "📄 Final Project Report":
    st.title("📄 Final Project Report")
    st.write("Prepare a short student project summary using the current molecule.")
    objective=st.text_area("Project objective","To explore molecular structures, visualize molecules in 2D and 3D, and calculate molecular descriptors.")
    observations=st.text_area("Observations","")
    conclusion=st.text_area("Conclusion","")
    if st.button("Generate project report"):
        report=f"""CHEMINFORMATICS VIRTUAL LABORATORY\nStudent: {student or 'Not entered'}\nRegistration: {regno or 'Not entered'}\nMolecule: {st.session_state.mol_name}\nSMILES: {st.session_state.smiles}\n\nOBJECTIVE\n{objective}\n\nOBSERVATIONS\n{observations}\n\nCONCLUSION\n{conclusion}\n"""
        st.download_button("Download project report (TXT)",report,"final_project_report.txt","text/plain")
        st.success("Project report prepared.")
