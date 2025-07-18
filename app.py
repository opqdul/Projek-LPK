koimport streamlit as st

# Kamus gugus fungsi
gugus_fungsi_kamus = {
    'COOH': 'Asam Karboksilat',
    'CHO': 'Aldehid',
    'CO': 'Keton',
    'OH': 'Alkohol',
    'NH2': 'Amina',
    'COO': 'Ester',
    'C=C': 'Alkena',
    'C≡C': 'Alkina'
}


# Kamus nama senyawa lengkap (tanpa strip)
kamus_nama_senyawa = {
    'CH3COOH': {'iupac': 'Asam etanoat', 'trivial': 'Asam asetat'},
    'HCOOH': {'iupac': 'Asam metanoat', 'trivial': 'Asam format'},
     'CH3OH': {'iupac': 'Metanol', 'trivial': 'Metil alkohol'},
    'CH3CH2OH': {'iupac': 'Etanol', 'trivial': 'Alkohol etil'},
    'CH3CHO': {'iupac': 'Etanal', 'trivial': 'Asetaldehida'},
    'CH3COCH3': {'iupac': 'Propanon', 'trivial': 'Aseton'},
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
    'CH4': {'iupac': 'Metana', 'trivial': '-'},
    'CH3CH3': {'iupac': 'Etana', 'trivial': '-'},
    'CH3CH2CH3': {'iupac': 'Propana', 'trivial': '-'},
    'CH3CH2CH2CH3': {'iupac': 'Butana', 'trivial': '-'},
    'CH3CH2CH2CH2CH3': {'iupac': 'Pentana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH3': {'iupac': 'Heksana', 'trivial': '-'}
    'CH4': {'iupac': 'Metana', 'trivial': '-'},
    'CH3CH3': {'iupac': 'Etana', 'trivial': '-'},
 'CH3CH2CH3': {'iupac': 'Propana', 'trivial': '-'},
    'CH3CH2CH2CH3': {'iupac': 'Butana', 'trivial': '-'},
    'CH3CH2CH2CH2CH3': {'iupac': 'Pentana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH3': {'iupac': 'Heksana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH3': {'iupac': 'Heptana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Oktana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonana', 'trivial': '-'},
    'CH3CH2CH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekana', 'trivial': '-'},
    'CH2=CH2': {'iupac': 'Etena', 'trivial': 'Etilena'},
    'CH2=CHCH3': {'iupac': 'Propena', 'trivial': 'Propilena'},
    'CH2=CHCH2CH3': {'iupac': 'Butena', 'trivial': '-'},
    'CH2=CHCH2CH2CH3': {'iupac': 'Pentena', 'trivial': '-'},
    'CH2=CHCH2CH2CH2CH3': {'iupac': 'Heksena', 'trivial': '-'},
    'CH2=CHCH2CH2CH2CH2CH3': {'iupac': 'Heptena', 'trivial': '-'},
    'CH2=CHCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktena', 'trivial': '-'},
    'CH2=CHCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonena', 'trivial': '-'},
    'CH2=CHCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekena', 'trivial': '-'}
 'CH≡CH': {'iupac': 'Etuna', 'trivial': 'Asetilena'},
    'CH≡CCH3': {'iupac': 'Propuna', 'trivial': '-'},
    'CH≡CCH2CH3': {'iupac': 'Butuna', 'trivial': '-'},
    'CH≡CCH2CH2CH3': {'iupac': 'Pentuna', 'trivial': '-'},
    'CH≡CCH2CH2CH2CH3': {'iupac': 'Heksuna', 'trivial': '-'},
    'CH≡CCH2CH2CH2CH2CH3': {'iupac': 'Heptuna', 'trivial': '-'},
    'CH≡CCH2CH2CH2CH2CH2CH3': {'iupac': 'Oktena', 'trivial': '-'},
    'CH≡CCH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Nonuna', 'trivial': '-'},
    'CH≡CCH2CH2CH2CH2CH2CH2CH2CH3': {'iupac': 'Dekuna', 'trivial': '-'},
}

# Fungsi identifikasi gugus fungsi
def identifikasi_gugus_fungsi(rumus):
    hasil = []
    for gugus, nama in gugus_fungsi_kamus.items():
        if gugus in rumus:
            hasil.append(nama)
    return hasil if hasil else ['Tidak teridentifikasi']

# Judul
st.title("🧪 Identifikasi Gugus Fungsi & Penamaan Senyawa Organik")

# Input
input_rumus = st.text_input("Masukkan rumus senyawa (contoh: CH3CH2OH atau CH3-CH2-OH):")

if input_rumus:
    # Normalisasi: hapus strip agar cocok dengan kamus
    rumus = input_rumus.replace("-", "").replace("=", "").replace("≡", "")

 hasil = identifikasi_gugus_fungsi(rumus)

    nama_iupac = "-"
    nama_trivial = "-"

    if rumus in kamus_nama_senyawa:
        nama_iupac = kamus_nama_senyawa[rumus]['iupac']
        nama_trivial = kamus_nama_senyawa[rumus]['trivial']
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
    st.write(f"*Rumus Diberikan:* {input_rumus}")
    st.write(f"*Rumus Distandarisasi:* {rumus}")
    st.write(f"*Gugus Fungsi Terdeteksi:* {', '.join(hasil)}")
    st.write(f"*Nama IUPAC:* {nama_iupac}")
    st.write(f"*Nama Trivial:* {nama_trivial}")
