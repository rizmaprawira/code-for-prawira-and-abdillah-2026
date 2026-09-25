# NCEP-R2 DJF 6-Experiment Moist LBM Matrix

Direktori ini berisi matriks komprehensif 6 eksperimen Linear Baroclinic Model (LBM) moist berbasis sirkulasi atmosfer musim dingin boreal (DJF) NCEP-R2 dan anomali SST regresi ERA5 terhadap indeks Niño-3.4.

Struktur direktori ini bersifat mandiri (*100% self-contained*), menggabungkan data basic state, data forcing (global & regional bersarang), skrip pre-processing, running, dan post-processing ke dalam satu wadah tanpa dependensi penamaan legasi (*legacy-free*).

---

## 📊 Matriks Eksperimen

| Nama Eksperimen | Basic State | Region 1 (C-E Pacific) | Region 2 (W Pacific) | Region 3 (Indian Ocean) | Deskripsi Fisika |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`CTRL`** | **P1** | P1 | P1 | P1 | Eksperimen kontrol referensi dengan basic state P1 (1981–2006) dan SST forcing P1 di seluruh domain. |
| **`EXP_B`** | **P2** | P1 | P1 | P1 | Mengisolasi efek perubahan basic state atmosfer: basic state diganti menjadi P2 (2007–2025) dengan mempertahankan SST forcing P1. |
| **`EXP_S-123`** | **P1** | P2 | P2 | P2 | Mengisolasi efek perubahan total SST forcing: forcing diganti menjadi P2 di ketiga wilayah dengan mempertahankan basic state P1. |
| **`EXP_B_S-123`** | **P2** | P2 | P2 | P2 | Konfigurasi P2 penuh: baik basic state maupun SST forcing di seluruh wilayah menggunakan periode P2. |
| **`EXP_B_S-1`** | **P2** | P2 | P1 | P1 | Menguji kontribusi perubahan SST forcing di Region 1 (Pasifik Tengah-Timur) di bawah latar belakang basic state P2; Region 2 & 3 tetap P1. |
| **`EXP_B_S-12`** | **P2** | P2 | P2 | P1 | Menguji kontribusi kumulatif perubahan SST forcing di Region 1 & Region 2 di bawah basic state P2; Region 3 tetap P1. |

### Definisi Koordinat Wilayah:
- **Region 1 (Pasifik Tengah–Timur)**: 160°E – 280°E, 20°S – 20°N
- **Region 2 (Pasifik Barat)**: 120°E – 160°E, 20°S – 20°N
- **Region 3 (Samudra Hindia)**: 30°E – 120°E, 20°S – 20°N

---

## 📁 Struktur Direktori

```text
ncep_r2_moist_experiment/
├── README_EXPERIMENTS.md                 # Dokumentasi ini
├── experiment_manifest.csv               # Manifest metadata ke-6 run
├── run_all_experiments.sh                # Skrip runner terpadu
├── combine_experiments_to_nc.sh          # Pipeline konversi CDO & kombinasi NetCDF
├── combine_experiments.py                # Skrip Python xarray penggabung NetCDF
├── gather_results.py                     # Skrip pengumpul plot & hasil
│
├── src/
│   ├── make_regional_forcing.py          # Generator forcing regional p2_s1 & p2_s12
│   └── setup_experiments.sh              # Skrip staging & stager mirror input lokal
│
├── inputs/                               # Master data induk (Self-contained)
│   ├── basic_state/                      # Basic state atmosfer & SST
│   │   ├── p1/                           # NCEP-R2 & ERA5 SST 1981-2006
│   │   ├── p2/                           # NCEP-R2 & ERA5 SST 2007-2025
│   │   └── shared/                       # grz.t21, gridx.t21, wgwin.t21
│   └── forcing/                          # Forcing SST regresi
│       ├── p1/                           # Global P1
│       ├── p2/                           # Global P2
│       ├── p2_s1/                        # P2 di Region 1, P1 di wilayah lain
│       └── p2_s12/                       # P2 di Region 1 & 2, P1 di wilayah lain
│
├── CTRL/                                 # Run directory 1
├── EXP_B/                                # Run directory 2
├── EXP_S-123/                            # Run directory 3
├── EXP_B_S-123/                          # Run directory 4
├── EXP_B_S-1/                            # Run directory 5
├── EXP_B_S-12/                           # Run directory 6
│
├── plot/                                 # Skrip visualisasi
│   └── plot_all_experiments.sh
└── results/                              # Hasil visualisasi t15 & t29
```

---

## 🚀 Panduan Eksekusi

### 1. Pre-Processing & Staging (Jika ingin re-stage)
```bash
# Generate ulang forcing regional p2_s1 & p2_s12 jika diperlukan
python3 src/make_regional_forcing.py

# Stage / refresh file input mirror ke seluruh 6 folder eksperimen
bash src/setup_experiments.sh
```

### 2. Validasi Dry-Run (Tanpa menjalankan model)
```bash
bash run_all_experiments.sh --dry-run
```

### 3. Eksekusi Model LBM
```bash
# Menjalankan 1 eksperimen spesifik:
bash run_all_experiments.sh --run CTRL

# Menjalankan seluruh 6 eksperimen berurutan:
bash run_all_experiments.sh --run-all

# Menjalankan seluruh 6 eksperimen dengan post-processing (gt2gr) tanpa plot otomatis:
bash run_all_experiments.sh --run-post-all
```

### 4. Post-Processing & Penggabungan Dataset NetCDF
```bash
# Menggabungkan seluruh output model ke 1 file NetCDF terpadu:
bash combine_experiments_to_nc.sh
# Output: combined_ncep_r2_moist_lbm.nc
```

### 5. Visualisasi & Plotting
```bash
# Menghasilkan plot t15 & t29 untuk psi850, psi500, z850, z500, dan psi_q:
bash plot/plot_all_experiments.sh

# Atau kumpulkan plot ke folder results/:
python3 gather_results.py
```
