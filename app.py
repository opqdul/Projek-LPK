import streamlit as st

# Kamus gugus fungsi
gugus_fungsi_kamus = {
    'COOH': 'Asam Karboksilat',
    'CHO': 'Aldehid',
    'CO': 'Keton',
    'OH': 'Alkohol',
    'NH2': 'Amina',
    'COO': 'Ester'    
}

ikatan = {
    'CHC': 'Alkuna',
    'CH2': 'Alkena',
    'CH':'Alkana'
}

# Kamus nama senyawa lengkap (tanpa strip)
kamus_nama_senyawa = {
    'HCOOCH3':{'iupac':'Metil metanoat','trivial':'',"golongan":"Ester","rumus_umum":"R-COO-R"},
    'HCOOC2H5': {'iupac': 'Etil metanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3COOCH3': {'iupac': 'Metil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3COOC2H5': {'iupac': 'Etil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3CH2COOCH3': {'iupac': 'Metil propanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3CH2COOC2H5': {'iupac': 'Etil propanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3CH2CH2COOCH3': {'iupac': 'Metil butanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3CH2CH2COOC2H5': {'iupac': 'Etil butanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3COOCH2CH2CH3': {'iupac': 'Propanil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3COOCH(CH3)2': {'iupac': 'Isopropil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'C6H5CH2COOCH3': {'iupac': 'Benzil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3CH2CH2COOCH2CH2CH3': {'iupac': 'Propanil butanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum': 'R-COO-R'},
    'CH3COOCH(CH2CH3)CH3': {'iupac': 'Oktil etanoat', 'trivial': '', 'golongan': 'Ester', 'rumus_umum':'R-COO-R'},
    'CH3NH2':{'iupac':'Metilamin','trivial':'','golongan':'Amina primer','rumus_umum':'R-NH2'},
    'CH3CH32NH': {'iupac': 'Dimetilamin', 'trivial': '', 'golongan': 'Amina sekunder', 'rumus_umum': 'R2NH'},
    'CH3CH3CH3N': {'iupac': 'Trimetilamin', 'trivial': '', 'golongan': 'Amina tersier', 'rumus_umum': 'R3N'},
    'CH3CH2NH2': {'iupac': 'Etilamin', 'trivial': '', 'golongan': 'Amina primer', 'rumus_umum': 'R-NH2'},
    'C6H5CH2NH2': {'iupac': 'Benzilamin', 'trivial': '', 'golongan': 'Amina primer', 'rumus_umum': 'ArCH2NH2'},
    'C6H5NH2': {'iupac': 'Anilin', 'trivial': '', 'golongan': 'Amina primer aromatik', 'rumus_umum': 'Ar-NH2'},
    'CH3OH': {'iupac': 'Metanol', 'trivial': 'Alkohol metil / spiritus', 'golongan': 'Alkohol primer', 'rumus_umum': 'R-OH'},
    'C2H5OH': {'iupac': 'Etanol', 'trivial': 'Alkohol etil / alkohol', 'golongan': 'Alkohol primer', 'rumus_umum': 'R-OH'},
    'CH3CH2CH2OH': {'iupac': '1-Propanol', 'trivial': 'n-Propanol', 'golongan': 'Alkohol primer', 'rumus_umum': 'R-OH'},
    'CH3CHOHCH3': {'iupac': '2-Propanol', 'trivial': 'Isopropanol / alkohol gosok', 'golongan': 'Alkohol sekunder', 'rumus_umum': 'R-OH'},
    'CH3CH2CH2CH2OH': {'iupac': '1-Butanol', 'trivial': 'n-Butanol', 'golongan': 'Alkohol primer', 'rumus_umum': 'R-OH'},
    'CH3CH2CHOHCH3': {'iupac': '2-Butanol', 'trivial': 'sec-Butanol', 'golongan': 'Alkohol sekunder', 'rumus_umum': 'R-OH'},
    '(CH3)2CHCH2OH': {'iupac': '2-Metil-1-propanol', 'trivial': 'Isobutanol', 'golongan': 'Alkohol primer', 'rumus_umum': 'R-OH'},
    '(CH3)3COH': {'iupac': '2-Metil-2-propanol', 'trivial': 't-Butanol', 'golongan': 'Alkohol tersier', 'rumus_umum': 'R-OH'},
    'C6H5CH2OH': {'iupac': 'Benzil alkohol', 'trivial': 'Benzil alkohol', 'golongan': 'Alkohol primer aromatik', 'rumus_umum': 'Ar-CH2OH'},
    'HOCH2CH2OH': {'iupac': 'Etana-1,2-diol', 'trivial': 'Glikol etilen', 'golongan': 'Alkohol diol', 'rumus_umum': 'HO-R-OH'},
    'HOCH2CH(OH)CH2OH': {'iupac': 'Propana-1,2,3-triol', 'trivial': 'Gliserol / gliserin', 'golongan': 'Alkohol triol', 'rumus_umum': 'R-(OH)3'},
    'CH3CH(OH)CH3': {'iupac': 'Propan-2-ol', 'trivial': 'Isopropil alkohol', 'golongan': 'Alkohol sekunder', 'rumus_umum':'R-OH'},
    'CH3COCH3': {'iupac': 'Propanon', 'trivial': 'Aseton', 'golongan': 'Keton', 'rumus_umum': 'R-CO-R'},
    'CH3COCH2CH3': {'iupac': 'Butanon', 'trivial': 'Metil etil keton', 'golongan': 'Keton', 'rumus_umum': 'R-CO-R'},
    'CH3COCH2CH2CH3': {'iupac': '2-Pentanon', 'trivial': 'Metil propil keton', 'golongan': 'Keton', 'rumus_umum': 'R-CO-R'},
    'CH3COCH2CH2CH2CH3': {'iupac': '2-Heksanon', 'trivial': 'Metil butil keton', 'golongan': 'Keton', 'rumus_umum': 'R-CO-R'},
    'C6H5COCH3': {'iupac': '1-Fenil-etanon', 'trivial': 'Asetofenon', 'golongan': 'Keton aromatik', 'rumus_umum': 'Ar-CO-R'},
    'C6H5COC6H5': {'iupac': '1,2-Difenil-etanon', 'trivial': 'Benzofenon', 'golongan': 'Keton aromatik', 'rumus_umum':'Ar-CO-Ar'},
    'HCHO': {'iupac': 'Metanal', 'trivial': 'Formaldehida', 'golongan': 'Aldehid', 'rumus_umum': 'R-CHO'},
    'CH3CHO': {'iupac': 'Etanal', 'trivial': 'Asetaldehida', 'golongan': 'Aldehid', 'rumus_umum': 'R-CHO'},
    'CH3CH2CHO': {'iupac': 'Propanal', 'trivial': 'Propionaldehida', 'golongan': 'Aldehid', 'rumus_umum': 'R-CHO'},
    'CH3(CH2)2CHO': {'iupac': 'Butanal', 'trivial': 'Butiraldehida', 'golongan': 'Aldehid', 'rumus_umum': 'R-CHO'},
    'CH3(CH2)3CHO': {'iupac': 'Pentanal', 'trivial': 'Valeraldehida', 'golongan': 'Aldehid', 'rumus_umum': 'R-CHO'},
    'C6H5CHO': {'iupac': 'Benzena karbaldehida', 'trivial': 'Benzaldehida', 'golongan': 'Aldehid aromatik', 'rumus_umum':'Ar-CHO'},
    'HCOOH': {'iupac': 'Asam metanoat', 'trivial': 'Asam format', 'golongan': 'Asam karboksilat', 'rumus_umum': 'R-COOH'},
    'CH3COOH': {'iupac': 'Asam etanoat', 'trivial': 'Asam asetat', 'golongan': 'Asam karboksilat', 'rumus_umum': 'R-COOH'},
    'CH3CH2COOH': {'iupac': 'Asam propanoat', 'trivial': 'Asam propionat', 'golongan': 'Asam karboksilat', 'rumus_umum': 'R-COOH'},
    'CH3(CH2)2COOH': {'iupac': 'Asam butanoat', 'trivial': 'Asam butirat', 'golongan': 'Asam karboksilat', 'rumus_umum': 'R-COOH'},
    'CH3(CH2)3COOH': {'iupac': 'Asam pentanoat', 'trivial': 'Asam valerianat', 'golongan': 'Asam karboksilat', 'rumus_umum': 'R-COOH'},
    'C6H5COOH': {'iupac': 'Asam benzoat', 'trivial': 'Asam benzoat', 'golongan': 'Asam karboksilat aromatik', 'rumus_umum':'Ar-COOH'},
    'CH4': {'iupac': 'Metana', 'trivial': '-','gambar':'metana.jpg'"golongan":,"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH3': {'iupac': 'Etana', 'trivial': '-','gambar':'etana.jpg'"golongan":,"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH3': {'iupac': 'Propana', 'trivial': '-','gambar':'propana.jpg'"golongan":,"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH3': {'iupac': 'Butana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH3': {'iupac': 'Pentana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH2CH3': {'iupac': 'Heksana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH2CH2CH3': {'iupac': 'Heptana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Oktana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekana', 'trivial': '-',"golongan":"Alkana","rumus_umum":"CnH2n+2"},
    'CH2CH2': {'iupac': 'Etena', 'trivial': 'Etilena',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH3': {'iupac': 'Propena', 'trivial': 'Propilena',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH3': {'iupac': 'Butena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH3': {'iupac': 'Pentena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH2CH3': {'iupac': 'Heksena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH2CH2CH3': {'iupac': 'Heptena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CH2CHCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekena', 'trivial': '-',"golongan":"Alkena","rumus_umum":"CnH2n"},
    'CHCH': {'iupac': 'Etuna', 'trivial': 'Asetilena',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH3': {'iupac': 'Propuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH3': {'iupac': 'Butuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH3': {'iupac': 'Pentuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH2CH3': {'iupac': 'Heksuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH2CH2CH3': {'iupac': 'Heptuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
    'CHCCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekuna', 'trivial': '-',"golongan":"Alkuna","rumus_umum":"CnH2n-2"},
}

def tentang():
    st.title("Isi tentang")

# Fungsi identifikasi gugus fungsi
def identifikasi_gugus_fungsi(rumus):
    hasil = []
    for gugus, nama in gugus_fungsi_kamus.items():
        if gugus in rumus:
            hasil.append(nama)
    return hasil if hasil else ['Tidak teridentifikasi']

def identifikasi_ikatan(rumus):
    for gugus, nama in ikatan.items():
        if rumus.find(gugus) == 0:
           return nama        
    return "Tidak teridentifikasi"


# Judul
def identifikasi():
    st.markdown(
        f"""
        <style>
        .blurimg {{
            background-image: url("https://raw.githubusercontent.com/RIVI44/LPK-KEDUA-/main/bg2.jpg");
            background-size: cover;
            position:fixed;
            inset:0;
            filter:blur(5px);
        }}        
        </style>
        <div class="blurimg"></div>
        """,
        unsafe_allow_html=True
    )
    
    st.title("🧪 Identifikasi Gugus Fungsi & Penamaan Senyawa Organik")

    # st.image("https://raw.githubusercontent.com/RIVI44/LPK-KEDUA-/main/WhatsApp%20Image%202025-07-19%20at%2013.17.34_bfbfabba.jpg", use_container_width=True)

    # Input
    input_rumus = st.text_input("Masukkan rumus senyawa (contoh: CH3CH2OH atau CH3-CH2-OH):")

    if input_rumus:
        # Normalisasi: hapus strip agar cocok dengan kamus
        rumus = input_rumus.replace("-", "").replace("=", "").replace("≡", "")

        hasil = identifikasi_gugus_fungsi(rumus)
        ikatan = identifikasi_ikatan(rumus)        

        nama_iupac = "-"
        nama_trivial = "-"
        gambar = None

        if rumus in kamus_nama_senyawa:
            nama_iupac = kamus_nama_senyawa[rumus]['iupac']
            nama_trivial = kamus_nama_senyawa[rumus]['trivial']
            if "gambar" in kamus_nama_senyawa[rumus]:
                gambar = kamus_nama_senyawa[rumus]['gambar']
        else:
            # Deteksi otomatis nama IUPAC sederhana
            if 'Asam Karboksilat' in hasil:
                nama_iupac = f"Asam {rumus.lower()}"
            elif 'Aldehid' in hasil:
                nama_iupac = f"{rumus.lower()} - al"
            elif 'Keton' in hasil:
                nama_iupac = f"{rumus.lower()} - on"
            elif 'Alkohol' in hasil:
                nama_iupac = f"{rumus.lower()} - ol"
            elif 'Amina' in hasil:
                nama_iupac = f"{rumus.lower()} - amina"

        # Output
        st.markdown("### 🔍 Hasil Identifikasi")
        
        if gambar:
            st.image(f"https://raw.githubusercontent.com/RIVI44/LPK-KEDUA-/main/{gambar}", width=250)
            
        with st.container(border=True):
            st.write(f"*Rumus Diberikan:* {input_rumus}")
            st.write(f"*Rumus Distandarisasi:* {rumus}")
            st.write(f"*Gugus Fungsi Terdeteksi:* {', '.join(hasil)}")
            st.write(f"*Jenis Hidrokarbon:* {ikatan}")
            st.write(f"*Nama IUPAC:* {nama_iupac}")
            st.write(f"*Nama Trivial:* {nama_trivial}")


option = st.sidebar.radio(
    "Menu:",
    ("Indentifikasi Gugus Fungsi", "Tentang Aplikasi")
)
if option == "Indentifikasi Gugus Fungsi":
    identifikasi();
elif option == "Tentang Aplikasi":
    tentang()

