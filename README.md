# ChemInformatics Virtual Laboratory

A Streamlit educational app using PubChem PUG REST for compound lookup and RDKit for local chemistry calculations and molecular rendering.

## Features
- Screenshot-aligned Home page, Theory, sidebar student information and learning objectives
- Example molecules table and developer/faculty information section
- PubChem search by compound name or CID
- SMILES validation
- 2D structure rendering with RDKit
- Interactive 3D model using 3Dmol.js in the browser
- Molecular descriptors and downloadable PDF report
- Example molecule comparisons and CSV graph export
- Theory, MCQ assessment, student information, theme colours, project report

## Run locally
Python 3.11 is a good starting point.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud
1. Create a new GitHub repository.
2. Upload `app.py`, `requirements.txt`, and `README.md` to the repository root.
3. Open https://share.streamlit.io/ and sign in with GitHub.
4. Choose **Create app**, select your repository and branch, and set the main file path to `app.py`.
5. Deploy.

## Important notes
- PubChem is accessed live through its public PUG REST API. The app searches records on demand; it does not download or mirror the entire PubChem database. Database size and available records change over time, and API rate limits/temporary outages may occur.
- 2D structures and molecular descriptors are calculated locally with RDKit.
- 3D coordinates are generated locally by RDKit for the selected molecule. 3Dmol.js is loaded from its public CDN in the browser.
- Some molecules may fail 3D conformer generation. Very large molecules may take longer.
- The descriptors are computational estimates, not experimental measurements or medical advice.
- This is an educational starter project. Test it before using it for coursework or public deployment.
