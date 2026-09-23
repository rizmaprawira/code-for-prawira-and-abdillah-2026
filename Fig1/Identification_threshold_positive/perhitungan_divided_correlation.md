# Perhitungan Divided Correlation

Dokumen ini merangkum workflow pada:

- `notebooks/divided_correlation_nino3/dcorr_nino3_workflow.ipynb`
- `notebooks/divided_correlation_nino4/dcorr_nino4_workflow.ipynb`
- `notebooks/divided_correlation_nino12/dcorr_nino12_workflow.ipynb`
- `notebooks/divided_correlation_nino34/dcorr_nino34_workflow.ipynb`

Keempat notebook memakai alur yang sama. Perbedaannya hanya pada indeks ENSO yang dibaca:

- Nino3
- Nino4
- Nino1+2
- Nino3.4

## Domain

Domain analisis yang dipakai di notebook adalah:

$$
\lambda \in [90.0^\circ E, 152.5^\circ E], \qquad
\phi \in [-12.5^\circ, 12.5^\circ]
$$

Domain ini kemudian dipotong lagi oleh mask daratan Natural Earth sehingga analisis sensitivitas utama difokuskan pada grid darat.

## Data Input

Input yang dipakai:

1. Data hujan bulanan MSWEP dari file NetCDF `mswep_monthly_combined.nc`.
2. Indeks ENSO bulanan dari CSV sesuai notebook:
   - `nino3.anom.csv`
   - `nino4.anom.csv`
   - `nina1.anom.csv` untuk notebook Nino1+2
   - `nino34.anom.csv`

Rentang waktu yang dibaca dari data bulanan adalah 1980-01 sampai 2025-12, lalu dibentuk menjadi DJF 1981 sampai DJF 2025.

## Output

Output yang dihasilkan notebook:

1. NetCDF hasil full-period correlation dan regression slope.
2. NetCDF hasil split-period maps:
   - `corr_p1`
   - `corr_p2`
   - `delta_corr = corr_p2 - corr_p1`
   - `slope_p1`
   - `slope_p2`
   - `delta_slope = slope_p2 - slope_p1`
3. CSV tabel split period, metrics, dan hotspot.
4. PNG peta baseline, split, atlas, summary, dan running diagnostics.
5. JSON manifest run.

## Algoritma Step by Step

### 1. Load data

Misalkan:

- $P(t, \lambda, \phi)$ = curah hujan bulanan MSWEP
- $N(t)$ = indeks ENSO bulanan

Data dibaca lalu diseragamkan koordinatnya. Untuk hujan, notebook memakai satuan mm/bulan.

### 2. Bentuk DJF seasonal mean

DJF untuk tahun $y$ didefinisikan sebagai:

$$
X_{\mathrm{DJF}}(y) = \frac{1}{3}
\left(
X_{\mathrm{Dec}(y-1)} +
X_{\mathrm{Jan}(y)} +
X_{\mathrm{Feb}(y)}
\right)
$$

Formula ini dipakai untuk hujan dan untuk indeks ENSO.

Hanya DJF year yang lengkap, yaitu yang punya tiga bulan penuh, yang dipertahankan.

### 3. Hitung klimatologi

Klimatologi dihitung pada periode 1991-2020.

Untuk hujan:

$$
\mu_P(\lambda, \phi) = \frac{1}{|C|} \sum_{k \in C} P_{\mathrm{DJF}}(k, \lambda, \phi)
$$

Untuk indeks ENSO:

$$
\mu_N = \frac{1}{|C|} \sum_{k \in C} N_{\mathrm{DJF}}(k)
$$

dengan $C$ adalah himpunan tahun 1991-2020.

### 4. Hitung anomali

Anomali hujan:

$$
P'(t, \lambda, \phi) = P_{\mathrm{DJF}}(t, \lambda, \phi) - \mu_P(\lambda, \phi)
$$

Anomali indeks ENSO:

$$
N'(t) = N_{\mathrm{DJF}}(t) - \mu_N
$$

Notebook tidak memakai standardisasi dan tidak melakukan detrending tambahan.

### 5. Bangun mask valid dan mask darat

Notebook membentuk:

$$
M_{\mathrm{valid}}(\lambda, \phi) = \bigwedge_{t \in T} \mathrm{isfinite}\left(P'(t, \lambda, \phi)\right)
$$

$$
M_{\mathrm{land}}(\lambda, \phi) = \text{NaturalEarthLand}(\lambda, \phi)
$$

Mask analisis darat:

$$
M_{\mathrm{analysis}}(\lambda, \phi) = M_{\mathrm{land}}(\lambda, \phi) \land M_{\mathrm{valid}}(\lambda, \phi)
$$

Mask ocean:

$$
M_{\mathrm{ocean}}(\lambda, \phi) = \neg M_{\mathrm{land}}(\lambda, \phi) \land M_{\mathrm{valid}}(\lambda, \phi)
$$

### 6. Bentuk split periode

Periode analisis dibagi menjadi dua bagian:

- $P_1 = [1981, \dots, 1981 + k - 1]$
- $P_2 = [1981 + k, \dots, 2025]$

Nilai $k$ disapu dari 10 sampai 35 tahun, sehingga label split menjadi:

$$
10v35,\ 11v34,\ \dots,\ 35v10
$$

Tujuannya adalah menguji sensitivitas korelasi terhadap pembagian periode.

### 7. Hitung korelasi full-period

Untuk setiap grid $(\lambda, \phi)$, korelasi Pearson antara anomali hujan dan anomali ENSO dihitung sebagai:

$$
r(\lambda, \phi) =
\frac{\operatorname{Cov}\left(P'(\lambda, \phi), N'\right)}
{\sigma_{P'(\lambda, \phi)} \, \sigma_{N'}}
$$

atau ekuivalen:

$$
r(\lambda, \phi) =
\frac{\sum_{k=1}^{n} \left(P'_k(\lambda, \phi) - \overline{P'}(\lambda, \phi)\right)\left(N'_k - \overline{N'}\right)}
{\sqrt{
\sum_{k=1}^{n} \left(P'_k(\lambda, \phi) - \overline{P'}(\lambda, \phi)\right)^2
\sum_{k=1}^{n} \left(N'_k - \overline{N'}\right)^2
}}
$$

Notebook juga menghitung slope regresi:

$$
\beta(\lambda, \phi) = \frac{\operatorname{Cov}\left(P'(\lambda, \phi), N'\right)}{\operatorname{Var}(N')}
$$

### 8. Hitung correlation untuk P1 dan P2

Untuk setiap split:

$$
r_{P_1}(\lambda, \phi), \beta_{P_1}(\lambda, \phi)
$$

dan

$$
r_{P_2}(\lambda, \phi), \beta_{P_2}(\lambda, \phi)
$$

Kemudian delta dihitung sebagai:

$$
\Delta r(\lambda, \phi) = r_{P_2}(\lambda, \phi) - r_{P_1}(\lambda, \phi)
$$

$$
\Delta \beta(\lambda, \phi) = \beta_{P_2}(\lambda, \phi) - \beta_{P_1}(\lambda, \phi)
$$

### 9. Hitung metrik sensitivitas

Notebook menghitung metrik utama pada grid darat.

#### a. Grid berubah pada ambang tertentu

Untuk ambang $\tau \in \{0.4, 0.6, 0.8\}$:

$$
n_{\tau} = \sum_{(\lambda,\phi) \in M_{\mathrm{analysis}}} \mathbf{1}\left(|\Delta r(\lambda, \phi)| > \tau\right)
$$

$$
f_{\tau} = \frac{n_{\tau}}{N_{\mathrm{land}}}
$$

dengan $N_{\mathrm{land}}$ = jumlah grid darat valid.

Notebook juga mencatat jumlah grid dengan delta positif:

$$
n^+_{\tau} = \sum_{(\lambda,\phi) \in M_{\mathrm{analysis}}} \mathbf{1}\left(\Delta r(\lambda, \phi) > \tau\right)
$$

#### b. Sign flip

Sebuah grid dihitung sebagai sign flip bila:

$$
\operatorname{sign}\left(r_{P_1}(\lambda, \phi)\right) \neq \operatorname{sign}\left(r_{P_2}(\lambda, \phi)\right)
$$

dan:

$$
|r_{P_1}(\lambda, \phi)| \ge 0.2, \qquad |r_{P_2}(\lambda, \phi)| \ge 0.2
$$

#### c. Strengthen / weaken

Penguatan:

$$
\Delta r(\lambda, \phi) > 0.2
$$

Pelemahan:

$$
\Delta r(\lambda, \phi) < -0.2
$$

#### d. Ringkasan split

Metrik ringkas yang dipakai untuk ranking split:

$$
\overline{|\Delta r|} = \frac{1}{N_{\mathrm{land}}} \sum_{(\lambda,\phi) \in M_{\mathrm{analysis}}} |\Delta r(\lambda, \phi)|
$$

Split kemudian diurutkan berdasarkan metrik ini.

### 10. Pilih hotspot

Notebook memilih $N=6$ hotspot dari split terbaik, yaitu split dengan $\overline{|\Delta r|}$ terbesar.

Hotspot adalah grid darat dengan $|\Delta r|$ terbesar, dengan jarak minimum antar hotspot sekitar 5 derajat.

Jika kandidat yang memenuhi jarak minimum tidak cukup, notebook mengisi sisa hotspot dari kandidat terkuat berikutnya.

### 11. Hitung running diagnostics pada hotspot

Untuk tiap hotspot, notebook menghitung statistik berjalan dengan jendela $w=15$ tahun.

Korelasi berjalan pada jendela pusat:

$$
r_t^{(w)} = \operatorname{corr}\left(
\{P'_{t-h}, \dots, P'_{t+h}\},
\{N'_{t-h}, \dots, N'_{t+h}\}
\right)
$$

dengan $h = \lfloor w/2 \rfloor$.

Notebook juga menghitung:

$$
\operatorname{Var}_t^{(w)}(P'), \qquad
\operatorname{Var}_t^{(w)}(N'), \qquad
\operatorname{Cov}_t^{(w)}(P', N')
$$

### 12. Simpan output

Hasil akhir disimpan sebagai:

- NetCDF full-period maps
- NetCDF split maps
- CSV split periods
- CSV split metrics
- CSV hotspots
- PNG baseline, splits, atlas, summaries, running stats
- JSON manifest run

## Ringkas Per Notebook

- `Nino3`: memakai `nino3.anom.csv`
- `Nino4`: memakai `nino4.anom.csv`
- `Nino1+2`: memakai `nina1.anom.csv`
- `Nino3.4`: memakai `nino34.anom.csv`

Secara alur hitung, keempatnya sama.
