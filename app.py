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
    'CH3COOH': {'iupac': 'Asam etanoat', 'trivial': 'Asam asetat'},
    'HCOOH': {'iupac': 'Asam metanoat', 'trivial': 'Asam format'},
    'CH3OH': {'iupac': 'Metanol', 'trivial': 'Alkohol kayu'},
    'CH3CH2OH': {'iupac': 'Etanol', 'trivial': 'Alkohol etil'},
    'HCOH': {'iupac': 'Metanal', 'trivial': 'Formaldehid'},
    'CH3CHO': {'iupac': 'Etanal', 'trivial': 'Asetaldehida'},
    'CH3NH2': {'iupac': 'Metilamina', 'trivial': '-'},
    'CH2CH2': {'iupac': 'Etena', 'trivial': 'Etilena'},
    'CH3CHCH2': {'iupac': 'Propena', 'trivial': 'Propilena'},
    'CH3CH2CHCH2': {'iupac': 'But-1-ena', 'trivial': '-'},
    'CH2CHCHCH2': {'iupac': 'Buta-1,3-diena', 'trivial': 'Butadiena'},
    'CH3CHCHCH3': {'iupac': 'But-2-ena', 'trivial': '-'},
    'CHCH': {'iupac': 'Etuna', 'trivial': 'Asetilena'},
    'CH3CCH': {'iupac': 'Propuna', 'trivial': '-'},
    'CH3CH2CCH': {'iupac': 'Butuna', 'trivial': '-'},
    'CH3CCCH3': {'iupac': 'Butuna', 'trivial': '-'},
    'CH4': {'iupac': 'Metana', 'trivial': '-','gambar':'metana.jpg'},
    'CH3CH3': {'iupac': 'Etana', 'trivial': '-','gambar':'etana.jpg'},
    'CH3CH2CH3': {'iupac': 'Propana', 'trivial': '-'},
    'CH3CH2CH2CH3': {'iupac': 'Butana', 'trivial': '-'},
    'CH3CH2CH2CH2CH3': {'iupac': 'Pentana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH3': {'iupac': 'Heksana', 'trivial': '-'},
    'CH3CH2CH3': {'iupac': 'Propana', 'trivial': '-'},
    'CH3CH2CH2CH3': {'iupac': 'Butana', 'trivial': '-'},
    'CH3CH2CH2CH2CH3': {'iupac': 'Pentana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH3': {'iupac': 'Heksana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH3': {'iupac': 'Heptana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Oktana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekana', 'trivial': '-'},
    'CH2CH2': {'iupac': 'Etena', 'trivial': 'Etilena'},
    'CH2CHCH3': {'iupac': 'Propena', 'trivial': 'Propilena'},
    'CH2CHCH2CH3': {'iupac': 'Butena', 'trivial': '-'},
    'CH2CHCH2CH2CH3': {'iupac': 'Pentena', 'trivial': '-'},
    'CH2CHCH2CH2CH2CH3': {'iupac': 'Heksena', 'trivial': '-'},
    'CH2CHCH2CH2CH2CH2CH3': {'iupac': 'Heptena', 'trivial': '-'},
    'CH2CHCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktena', 'trivial': '-'},
    'CH2CHCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonena', 'trivial': '-'},
    'CH2CHCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekena', 'trivial': '-'},
    'CHCH': {'iupac': 'Etuna', 'trivial': 'Asetilena'},
    'CHCCH3': {'iupac': 'Propuna', 'trivial': '-'},
    'CHCCH2CH3': {'iupac': 'Butuna', 'trivial': '-'},
    'CHCCH2CH2CH3': {'iupac': 'Pentuna', 'trivial': '-'},
    'CHCCH2CH2CH2CH3': {'iupac': 'Heksuna', 'trivial': '-'},
    'CHCCH2CH2CH2CH2CH3': {'iupac': 'Heptuna', 'trivial': '-'},
    'CHCCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktena', 'trivial': '-'},
    'CHCCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonuna', 'trivial': '-'},
    'CHCCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekuna', 'trivial': '-'},
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
       if gugus in rumus:
           return nama
    return "Tidak teridentifikasi"


# Judul
def identifikasi():
    st.title("🧪 Identifikasi Gugus Fungsi & Penamaan Senyawa Organik")
    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("https://raw.githubusercontent.com/RIVI44/LPK-KEDUA-/main/bg2.jpg");
            background-size: cover;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

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
            
        st.markdown(f":blue-background[*Rumus Diberikan:* {input_rumus}]")
        st.markdown(f":white-background[*Rumus Distandarisasi:* {rumus}]")
        st.markdown(f":white-background[*Gugus Fungsi Terdeteksi:* {', '.join(hasil)}]")
        st.markdown(f":white-background[*Jenis Hidrokarbon:* {ikatan}]")
        st.markdown(f":white-background[*Nama IUPAC:* {nama_iupac}]")
        st.markdown(f":white-background[*Nama Trivial:* {nama_trivial}]")


option = st.sidebar.radio(
    "Menu:",
    ("Indentifikasi Gugus Fungsi", "Tentang Aplikasi")
)
if option == "Indentifikasi Gugus Fungsi":
    identifikasi();
elif option == "Tentang Aplikasi":
    tentang()

