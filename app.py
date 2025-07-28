import streamlit as st
import pandas as pd
import re
from typing import Dict, List, Tuple, Optional
import base64
from io import BytesIO

# Try to import chemistry libraries, fall back gracefully
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    import pubchempy as pcp
    PUBCHEM_AVAILABLE = True
except ImportError:
    PUBCHEM_AVAILABLE = False

# Configure page
st.set_page_config(
    page_title="ChemID Pro",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better UI
st.markdown("""
<style>
    .main {
        padding-top: 2rem;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        padding-left: 20px;
        padding-right: 20px;
    }
    .result-container {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 20px;
        border-radius: 10px;
        color: white;
        margin: 10px 0;
    }
    .property-card {
        background: rgba(255, 255, 255, 0.1);
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
        backdrop-filter: blur(10px);
    }
    .search-container {
        background: #f8f9fa;
        padding: 20px;
        border-radius: 10px;
        margin: 10px 0;
    }
    .stAlert > div {
        padding: 1rem;
        border-radius: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

class CompoundDatabase:
    """Efficient compound database with optimized search capabilities"""
    
    def __init__(self):
        self.compounds = self._load_compounds()
        self.name_index = self._build_name_index()
        
    def _load_compounds(self) -> Dict:
        """Load compound data efficiently"""
        return {
            # Alkanes
            'CH4': {'iupac': 'Metana', 'trivial': 'Metana', 'smiles': 'C', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            'C2H6': {'iupac': 'Etana', 'trivial': 'Etana', 'smiles': 'CC', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            'C3H8': {'iupac': 'Propana', 'trivial': 'Propana', 'smiles': 'CCC', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            'C4H10': {'iupac': 'Butana', 'trivial': 'Butana', 'smiles': 'CCCC', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            'C5H12': {'iupac': 'Pentana', 'trivial': 'Pentana', 'smiles': 'CCCCC', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            'C6H14': {'iupac': 'Heksana', 'trivial': 'Heksana', 'smiles': 'CCCCCC', 'group': 'Alkana', 'formula_pattern': 'CnH2n+2'},
            
            # Alkenes
            'C2H4': {'iupac': 'Etena', 'trivial': 'Etilena', 'smiles': 'C=C', 'group': 'Alkena', 'formula_pattern': 'CnH2n'},
            'C3H6': {'iupac': 'Propena', 'trivial': 'Propilena', 'smiles': 'C=CC', 'group': 'Alkena', 'formula_pattern': 'CnH2n'},
            'C4H8': {'iupac': 'Butena', 'trivial': 'Butena', 'smiles': 'C=CCC', 'group': 'Alkena', 'formula_pattern': 'CnH2n'},
            
            # Alkynes
            'C2H2': {'iupac': 'Etuna', 'trivial': 'Asetilena', 'smiles': 'C#C', 'group': 'Alkuna', 'formula_pattern': 'CnH2n-2'},
            'C3H4': {'iupac': 'Propuna', 'trivial': 'Propuna', 'smiles': 'C#CC', 'group': 'Alkuna', 'formula_pattern': 'CnH2n-2'},
            
            # Alcohols
            'CH4O': {'iupac': 'Metanol', 'trivial': 'Alkohol metil', 'smiles': 'CO', 'group': 'Alkohol', 'formula_pattern': 'R-OH'},
            'C2H6O': {'iupac': 'Etanol', 'trivial': 'Alkohol etil', 'smiles': 'CCO', 'group': 'Alkohol', 'formula_pattern': 'R-OH'},
            'C3H8O': {'iupac': '1-Propanol', 'trivial': 'n-Propanol', 'smiles': 'CCCO', 'group': 'Alkohol', 'formula_pattern': 'R-OH'},
            'C3H8O_iso': {'iupac': '2-Propanol', 'trivial': 'Isopropanol', 'smiles': 'CC(C)O', 'group': 'Alkohol', 'formula_pattern': 'R-OH'},
            'C4H10O': {'iupac': '1-Butanol', 'trivial': 'n-Butanol', 'smiles': 'CCCCO', 'group': 'Alkohol', 'formula_pattern': 'R-OH'},
            
            # Aldehydes
            'CH2O': {'iupac': 'Metanal', 'trivial': 'Formaldehida', 'smiles': 'C=O', 'group': 'Aldehid', 'formula_pattern': 'R-CHO'},
            'C2H4O': {'iupac': 'Etanal', 'trivial': 'Asetaldehida', 'smiles': 'CC=O', 'group': 'Aldehid', 'formula_pattern': 'R-CHO'},
            'C3H6O_ald': {'iupac': 'Propanal', 'trivial': 'Propionaldehida', 'smiles': 'CCC=O', 'group': 'Aldehid', 'formula_pattern': 'R-CHO'},
            'C4H8O_ald': {'iupac': 'Butanal', 'trivial': 'Butiraldehida', 'smiles': 'CCCC=O', 'group': 'Aldehid', 'formula_pattern': 'R-CHO'},
            
            # Ketones
            'C3H6O': {'iupac': 'Propanon', 'trivial': 'Aseton', 'smiles': 'CC(=O)C', 'group': 'Keton', 'formula_pattern': 'R-CO-R'},
            'C4H8O_ket': {'iupac': 'Butanon', 'trivial': 'Metil etil keton', 'smiles': 'CCC(=O)C', 'group': 'Keton', 'formula_pattern': 'R-CO-R'},
            'C5H10O': {'iupac': '2-Pentanon', 'trivial': 'Metil propil keton', 'smiles': 'CCCC(=O)C', 'group': 'Keton', 'formula_pattern': 'R-CO-R'},
            
            # Carboxylic acids
            'CH2O2': {'iupac': 'Asam metanoat', 'trivial': 'Asam format', 'smiles': 'C(=O)O', 'group': 'Asam Karboksilat', 'formula_pattern': 'R-COOH'},
            'C2H4O2': {'iupac': 'Asam etanoat', 'trivial': 'Asam asetat', 'smiles': 'CC(=O)O', 'group': 'Asam Karboksilat', 'formula_pattern': 'R-COOH'},
            'C3H6O2': {'iupac': 'Asam propanoat', 'trivial': 'Asam propionat', 'smiles': 'CCC(=O)O', 'group': 'Asam Karboksilat', 'formula_pattern': 'R-COOH'},
            'C4H8O2': {'iupac': 'Asam butanoat', 'trivial': 'Asam butirat', 'smiles': 'CCCC(=O)O', 'group': 'Asam Karboksilat', 'formula_pattern': 'R-COOH'},
            
            # Esters
            'C3H6O2_est': {'iupac': 'Metil etanoat', 'trivial': 'Metil asetat', 'smiles': 'CC(=O)OC', 'group': 'Ester', 'formula_pattern': 'R-COO-R'},
            'C4H8O2_est': {'iupac': 'Etil etanoat', 'trivial': 'Etil asetat', 'smiles': 'CC(=O)OCC', 'group': 'Ester', 'formula_pattern': 'R-COO-R'},
            'C4H8O2_est2': {'iupac': 'Metil propanoat', 'trivial': 'Metil propionat', 'smiles': 'CCC(=O)OC', 'group': 'Ester', 'formula_pattern': 'R-COO-R'},
            
            # Amines
            'CH5N': {'iupac': 'Metilamin', 'trivial': 'Metilamin', 'smiles': 'CN', 'group': 'Amina primer', 'formula_pattern': 'R-NH2'},
            'C2H7N': {'iupac': 'Etilamin', 'trivial': 'Etilamin', 'smiles': 'CCN', 'group': 'Amina primer', 'formula_pattern': 'R-NH2'},
            'C3H9N': {'iupac': 'Propilamin', 'trivial': 'Propilamin', 'smiles': 'CCCN', 'group': 'Amina primer', 'formula_pattern': 'R-NH2'},
            
            # Aromatic compounds
            'C6H6': {'iupac': 'Benzena', 'trivial': 'Benzena', 'smiles': 'c1ccccc1', 'group': 'Aromatik', 'formula_pattern': 'Ar-H'},
            'C7H8': {'iupac': 'Toluena', 'trivial': 'Metilbenzena', 'smiles': 'Cc1ccccc1', 'group': 'Aromatik', 'formula_pattern': 'Ar-R'},
            'C7H6O': {'iupac': 'Benzaldehida', 'trivial': 'Benzaldehida', 'smiles': 'O=Cc1ccccc1', 'group': 'Aldehid aromatik', 'formula_pattern': 'Ar-CHO'},
        }
    
    def _build_name_index(self) -> Dict:
        """Build index for fast name-based search"""
        index = {}
        for formula, data in self.compounds.items():
            # Index by IUPAC name
            index[data['iupac'].lower()] = formula
            # Index by trivial name
            if data['trivial'] and data['trivial'] != '-' and data['trivial'] != data['iupac']:
                index[data['trivial'].lower()] = formula
        return index
    
    def search_by_formula(self, formula: str) -> Optional[Dict]:
        """Search compound by molecular formula"""
        formula = self._normalize_formula(formula)
        return self.compounds.get(formula)
    
    def search_by_name(self, name: str) -> Optional[Tuple[str, Dict]]:
        """Search compound by name (IUPAC or trivial)"""
        name_lower = name.strip().lower()
        formula = self.name_index.get(name_lower)
        if formula:
            return formula, self.compounds[formula]
        return None
    
    def _normalize_formula(self, formula: str) -> str:
        """Normalize molecular formula"""
        # Remove spaces and common separators
        formula = re.sub(r'[\s\-=≡]', '', formula)
        return formula
    
    def get_functional_groups(self, smiles: str) -> List[str]:
        """Identify functional groups from SMILES"""
        groups = []
        
        if not RDKIT_AVAILABLE:
            # Fallback pattern matching
            if 'O' in smiles and '=' not in smiles and '#' not in smiles and '(=O)' not in smiles:
                groups.append('Alkohol')
            if 'C=O' in smiles and 'OC' not in smiles:
                groups.append('Aldehid/Keton')
            if 'C(=O)O' in smiles:
                groups.append('Asam Karboksilat')
            if 'C(=O)O' in smiles and 'OC' in smiles:
                groups.append('Ester')
            if 'N' in smiles and '=' not in smiles:
                groups.append('Amina')
            if 'c1ccccc1' in smiles:
                groups.append('Aromatik')
            return groups or ['Hidrokarbon']
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return ['Tidak dikenal']
            
            # Check for functional groups using SMARTS patterns
            patterns = {
                'Alkohol': '[OH]',
                'Aldehid': '[CX3H1](=O)[#6]',
                'Keton': '[#6][CX3](=O)[#6]',
                'Asam Karboksilat': '[CX3](=O)[OX2H1]',
                'Ester': '[#6][CX3](=O)[OX2H0][#6]',
                'Amina': '[NX3;H2,H1;!$(NC=O)]',
                'Aromatik': 'c1ccccc1'
            }
            
            for group_name, pattern in patterns.items():
                if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                    groups.append(group_name)
            
            return groups or ['Hidrokarbon']
        except:
            return ['Tidak dikenal']

class MoleculeVisualizer:
    """Handle molecule visualization using RDKit or fallback"""
    
    @staticmethod
    def draw_molecule(smiles: str, size: Tuple[int, int] = (300, 300)) -> str:
        """Draw molecule and return as base64 string with fallback"""
        if not RDKIT_AVAILABLE:
            return None
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            img = Draw.MolToImage(mol, size=size)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            return f"data:image/png;base64,{img_str}"
        except Exception as e:
            # Debug: show error
            st.error(f"Error generating molecule image: {str(e)}")
            return None
    
    @staticmethod
    def generate_ascii_structure(formula: str) -> str:
        """Generate ASCII representation when RDKit is not available"""
        ascii_structures = {
            'CH4': """
        H
        |
    H-C-H
        |
        H
            """,
            'C2H5OH': """
        H   H
        |   |
    H-C-C-O-H
        |   |
        H   H
            """,
            'CH3OH': """
        H
        |
    H-C-O-H
        |
        H
            """,
            'CH3COOH': """
        H     O
        |     ||
    H-C-C-O-H
        |
        H
            """,
            'CH3CH3': """
        H   H
        |   |
    H-C-C-H
        |   |
        H   H
            """,
            'CH3CH2CH3': """
        H   H   H
        |   |   |
    H-C-C-C-H
        |   |   |
        H   H   H
            """
        }
        
        return ascii_structures.get(formula, f"Chemical Formula: {formula}")
    
    @staticmethod
    def calculate_molecular_weight(formula: str) -> float:
        """Simple molecular weight calculation"""
        elements = {
            'C': 12.01, 'H': 1.008, 'O': 15.999, 'N': 14.007, 
            'S': 32.06, 'P': 30.97, 'Cl': 35.45, 'Br': 79.90
        }
        
        total_weight = 0
        i = 0
        while i < len(formula):
            if formula[i].isupper():
                element = formula[i]
                i += 1
                # Check for lowercase letter (two-letter element)
                if i < len(formula) and formula[i].islower():
                    element += formula[i]
                    i += 1
                
                # Get the count
                count_str = ""
                while i < len(formula) and formula[i].isdigit():
                    count_str += formula[i]
                    i += 1
                
                count = int(count_str) if count_str else 1
                
                if element in elements:
                    total_weight += elements[element] * count
            else:
                i += 1
        
        return round(total_weight, 2)

class ChemicalAnalyzer:
    """Main chemical analysis engine"""
    
    def __init__(self):
        self.db = CompoundDatabase()
        self.visualizer = MoleculeVisualizer()
    
    def analyze_compound(self, query: str, search_type: str) -> Dict:
        """Analyze compound and return comprehensive results"""
        result = {
            'found': False,
            'query': query,
            'search_type': search_type,
            'compound_data': None,
            'formula': None,
            'functional_groups': [],
            'properties': {},
            'visualization': None
        }
        
        if search_type == "formula":
            compound_data = self.db.search_by_formula(query)
            if compound_data:
                result['found'] = True
                result['formula'] = query
                result['compound_data'] = compound_data
        else:
            search_result = self.db.search_by_name(query)
            if search_result:
                formula, compound_data = search_result
                result['found'] = True
                result['formula'] = formula
                result['compound_data'] = compound_data
        
        if result['found']:
            smiles = result['compound_data']['smiles']
            result['functional_groups'] = self.db.get_functional_groups(smiles)
            result['visualization'] = self.visualizer.draw_molecule(smiles)
            
            # Get additional properties
            if RDKIT_AVAILABLE:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        result['properties'] = {
                            'Berat Molekul': f"{Descriptors.MolWt(mol):.2f} g/mol",
                            'Jumlah Atom': mol.GetNumAtoms(),
                            'Jumlah Ikatan': mol.GetNumBonds(),
                            'LogP': f"{Descriptors.MolLogP(mol):.2f}"
                        }
                except:
                    pass
        
        return result

# Initialize the analyzer
@st.cache_resource
def get_analyzer():
    return ChemicalAnalyzer()

def main():
    """Main application function with modern UI"""

    # Header
    st.markdown("""
    <div style='text-align: center; padding: 2rem 0;'>
        <h1 style='color: #2E86AB; margin-bottom: 0.5rem;'>ChemID Pro</h1>
        <p style='color: #666; font-size: 1.2rem;'>Identifikasi Senyawa Kimia dengan Teknologi Modern</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Create tabs for different features
    tab1, tab2, tab3 = st.tabs(["Identifikasi Senyawa", "Database Browser", "Tentang"])
    
    with tab1:
        compound_identifier()
    
    with tab2:
        database_browser()
    
    with tab3:
        about_page()

def compound_identifier():
    """Modern compound identification interface"""
    
    analyzer = get_analyzer()
    
    st.markdown('<div class="search-container">', unsafe_allow_html=True)
    
    # Search options
    col1, col2 = st.columns([1, 2])
    
    with col1:
        search_type = st.selectbox(
            "Pilih metode pencarian:",
            ["formula", "name"],
            format_func=lambda x: "Rumus Molekul" if x == "formula" else "Nama Senyawa"
        )
    
    with col2:
        if search_type == "formula":
            query = st.text_input(
                "Masukkan rumus molekul:",
                placeholder="Contoh: CH4, C2H6, C2H5OH",
                help="Masukkan rumus molekul standar"
            )
        else:
            query = st.text_input(
                "Masukkan nama senyawa:",
                placeholder="Contoh: metana, etanol, asam asetat",
                help="Masukkan nama IUPAC atau nama trivial"
            )
    
    st.markdown('</div>', unsafe_allow_html=True)
    
    if query:
        with st.spinner("Menganalisis senyawa..."):
            result = analyzer.analyze_compound(query, search_type)
        
        if result['found']:
            display_compound_results(result)
        else:
            st.error("❌ Senyawa tidak ditemukan dalam database")
            st.info("Coba gunakan format yang berbeda atau periksa ejaan")

def display_compound_results(result: Dict):
    """Display compound analysis results with modern UI"""
    
    data = result['compound_data']
    formula = result['formula']
    
    st.markdown('<div class="result-container">', unsafe_allow_html=True)
    
    # Header with compound info
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown(f"""
        ### {data['iupac']}
        **Rumus Molekul:** `{formula}`  
        **Nama Trivial:** {data['trivial']}  
        **Golongan:** {data['group']}  
        **Pola Rumus:** {data['formula_pattern']}
        """)
    
    with col2:
        if result['visualization']:
            st.markdown("**Struktur 2D:**")
            st.markdown(f'<img src="{result["visualization"]}" width="200">', unsafe_allow_html=True)
        elif not RDKIT_AVAILABLE:
            st.markdown("**Struktur ASCII:**")
            analyzer = get_analyzer()
            ascii_structure = analyzer.visualizer.generate_ascii_structure(formula)
            st.code(ascii_structure, language=None)
            st.caption("Install RDKit untuk visualisasi 2D yang lebih baik")
        else:
            st.info("Struktur tidak dapat ditampilkan")
    
    st.markdown('</div>', unsafe_allow_html=True)
    
    # Properties section
    if result['functional_groups'] or result['properties']:
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown('<div class="property-card">', unsafe_allow_html=True)
            st.markdown("### Gugus Fungsi")
            for group in result['functional_groups']:
                st.markdown(f"- {group}")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col2:
            if result['properties']:
                st.markdown('<div class="property-card">', unsafe_allow_html=True)
                st.markdown("### Sifat Molekul")
                for prop, value in result['properties'].items():
                    st.markdown(f"**{prop}:** {value}")
                st.markdown('</div>', unsafe_allow_html=True)
            elif not RDKIT_AVAILABLE:
                st.markdown('<div class="property-card">', unsafe_allow_html=True)
                st.markdown("### Sifat Molekul")
                analyzer = get_analyzer()
                mw = analyzer.visualizer.calculate_molecular_weight(formula)
                st.markdown(f"**Berat Molekul:** {mw} g/mol")
                st.caption("Install RDKit untuk analisis properti lengkap")
                st.markdown('</div>', unsafe_allow_html=True)

def database_browser():
    """Browse available compounds in database"""
    
    analyzer = get_analyzer()
    
    st.markdown("### 📚 Jelajahi Database Senyawa")
    
    # Filter options
    col1, col2 = st.columns(2)
    
    with col1:
        group_filter = st.selectbox(
            "Filter berdasarkan golongan:",
            ["Semua"] + list(set(data['group'] for data in analyzer.db.compounds.values()))
        )
    
    with col2:
        search_term = st.text_input("Cari senyawa:", placeholder="Ketik nama atau rumus...")
    
    # Display compounds
    compounds_to_show = []
    
    for formula, data in analyzer.db.compounds.items():
        # Apply filters
        if group_filter != "Semua" and data['group'] != group_filter:
            continue
        
        if search_term and search_term.lower() not in formula.lower() and \
           search_term.lower() not in data['iupac'].lower() and \
           search_term.lower() not in data['trivial'].lower():
            continue
        
        compounds_to_show.append({
            'Rumus': formula,
            'Nama IUPAC': data['iupac'],
            'Nama Trivial': data['trivial'],
            'Golongan': data['group']
        })
    
    if compounds_to_show:
        df = pd.DataFrame(compounds_to_show)
        
        # Make the table interactive
        event = st.dataframe(
            df,
            use_container_width=True,
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Show details for selected compound
        if event.selection.rows:
            selected_idx = event.selection.rows[0]
            selected_formula = compounds_to_show[selected_idx]['Rumus']
            
            st.markdown("---")
            st.markdown("### Detail Senyawa Terpilih")
            
            result = analyzer.analyze_compound(selected_formula, "formula")
            if result['found']:
                display_compound_results(result)
    else:
        st.info("Tidak ada senyawa yang cocok dengan filter yang dipilih.")

def about_page():
    """About page with modern design"""
    
    st.markdown("""
    ### Tentang ChemID Pro
    
    **ChemID Pro** adalah aplikasi modern untuk identifikasi dan analisis senyawa kimia organik. 
    Aplikasi ini dirancang untuk membantu pelajar, mahasiswa, dan peneliti dalam:
    
    #### Fitur Utama:
    - **Pencarian Dual-Mode**: Cari berdasarkan rumus molekul atau nama senyawa
    - **Visualisasi 2D**: Struktur molekul interaktif menggunakan RDKit
    - **Analisis Properti**: Berat molekul, LogP, dan karakteristik lainnya
    - **Deteksi Gugus Fungsi**: Identifikasi otomatis gugus fungsi
    - **Database Terintegrasi**: Akses ke database senyawa yang komprehensif
    
    #### Teknologi:
    - **Streamlit**: Framework aplikasi web modern
    - **RDKit**: Library cheminformatics terdepan
    - **PubChemPy**: Akses ke database PubChem
    - **Pandas**: Manipulasi dan analisis data
    
    #### Tujuan Pendidikan:
    Aplikasi ini dikembangkan untuk mendukung pembelajaran kimia organik dengan:
    - Interface yang intuitif dan user-friendly
    - Visualisasi yang membantu pemahaman struktur molekul
    - Data yang akurat dan terpercaya
    - Akses mudah ke informasi senyawa kimia

    """, unsafe_allow_html=True)

# Run the app
if __name__ == "__main__":
    main()
