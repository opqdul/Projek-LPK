# filename: app.py
import streamlit as st

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

# Kamus nama IUPAC dan trivial untuk senyawa umum (bisa dikembangkan)
kamus_nama_senyawa = {
    'CH3COOH': {
        'iupac': 'Asam etanoat',
        'trivial': 'Asam asetat'
    },
    'HCOOH': {
        'iupac': 'Asam metanoat',
        'trivial': 'Asam format'
    },
    'CH3OH': {
        'iupac': 'Metanol',
        'trivial': 'Alkohol kayu'
    },
    'CH3CH2OH': {
        'iupac': 'Etanol',
        'trivial': 'Alkohol etil'
    },
    'CH3CHO': {
        'iupac': 'Etanal',
        'trivial': 'Asetaldehida'
    },
    'CH3COCH3': {
        'iupac': 'Propanon',
        'trivial': 'Aseton'
    },
    'CH3NH2': {
        'iupac': 'Metilamina',
        'trivial': '-'
    },
    'CH2=CH2': {
        'iupac': 'Etena',
        'trivial': 'Etilena'
    },
    'CH≡CH': {
        'iupac': 'Etuna',
        'trivial': 'Asetilena'
    }
}

# Fungsi identifikasi gugus fungsi
def identifikasi_gugus_fungsi(rumus):
    hasil = []
    for gugus, nama in gugus_fungsi_kamus.items():
        if gugus in rumus:
            hasil.append(nama)
    return hasil if hasil else ['Tidak teridentifikasi']

# Judul halaman
st.title("🧪 Identifikasi Gugus Fungsi & Penamaan Senyawa Organik")

# Input user
rumus = st.text_input("Masukkan rumus senyawa (contoh: CH3COOH):")

if rumus:
    hasil = identifikasi_gugus_fungsi(rumus)

    # Nama senyawa default
    nama_iupac = "-"
    nama_trivial = "-"

    if rumus in kamus_nama_senyawa:
        nama_iupac = kamus_nama_senyawa[rumus]['iupac']
        nama_trivial = kamus_nama_senyawa[rumus]['trivial']
    else:
        # Jika tidak ada di kamus, buat nama IUPAC sangat sederhana
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

    # Output hasil
    st.markdown("### 🔍 Hasil Identifikasi")
    st.write(f"*Rumus Senyawa:* {rumus}")
    st.write(f"*Gugus Fungsi Terdeteksi:* {', '.join(hasil)}")
    st.write(f"*Nama IUPAC:* {nama_iupac}")
    st.write(f"*Nama Trivial:* {nama_trivial}")
