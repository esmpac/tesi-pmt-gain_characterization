# tesi-pmt-gain_characterization




This project implements a full analysis pipeline for the characterization of Photomultiplier Tubes (PMTs) used in high-energy physics experiments.

The goal is to determine the PMT gain as a function of the applied high voltage (HV), through a complete calibration workflow including pedestal subtraction, single-photoelectron (SPE) calibration, and gain curve reconstruction.

The analysis is based on ROOT and follows standard detector calibration techniques used in astroparticle physics.

---

## Data

The dataset consists of experimental PMT acquisition runs under different operating conditions:

- **Pedestal runs (dark conditions)**  
  Used to estimate the electronic baseline (ADC offset) for each channel.

- **Single-photoelectron (SPE) runs**  
  Low-light data used to calibrate the charge response in terms of photoelectrons.

- **Gain curve runs (high illumination)**  
  Data acquired at different high voltages (900–1500 V) used to reconstruct the gain curve.

All raw data are stored in CSV format and processed using ROOT-based histograms.

---

## Requirements

- ROOT (v6+)
- numpy
- pandas
- pathlib
- IPython (Jupyter/Colab visualization)

---

## Analysis Pipeline

The analysis is structured into three main stages:

---

### 1. Pedestal Estimation

- Extraction of the electronic baseline per channel  
- Gaussian fit of ADC distributions under dark conditions  
- Estimation of the pedestal mean:

$$
\mu_{\text{ped}}
$$

---

### 2. Single-Photoelectron (SPE) Calibration

- Pedestal subtraction from low-light data  
- Modeling of the charge spectrum using a double-Gaussian function  
- Extraction of the SPE mean charge per channel  
- Conversion from ADC units to physical charge (pC):

$$
Q_{\text{SPE}} = \mu_{\text{SPE}} \cdot C_{\text{ADC}}
$$

where \( C_{\text{ADC}} \) is the ADC-to-charge calibration constant.

- Conversion to number of photoelectrons:

$$
N_{pe} = \frac{Q_{\text{SPE}}}{e}
$$

---

### 3. Gain Curve Reconstruction

- Analysis of charge distributions as a function of high voltage (HV)  
- Gaussian fit of high-intensity charge spectra  
- Extraction of collected charge \( Q(V) \) at each voltage  
- Conversion to gain using SPE calibration:

$$
G(V) = \frac{Q(V)}{N_{pe} \cdot e}
$$

---

### Empirical Gain Model

The PMT gain is expected to follow a power-law dependence on the applied voltage:

$$
G(V) = A \cdot V^{b}
$$

where:
- \( A \) is a normalization constant  
- \( b \) is the gain exponent characteristic of the PMT structure

---

## Physical Quantities Extracted

The pipeline allows estimation of:

- Electronic pedestal (ADC baseline)  
- Single-photoelectron charge response  
- Number of photoelectrons per channel  
- PMT gain as a function of HV  
- Power-law exponent \( b \) of gain curves  
- Characteristic operating voltage per channel  

---

## Repository Structure

code/ → main ROOT-based analysis pipeline
data/ → raw experimental datasets (excluded via .gitignore)
plots/ → generated figures (gain curves, SPE fits, pedestal distributions)
root/ → optional ROOT outputs



---

## Notes

- Stable PMT response is assumed across acquisition runs  
- Accurate pedestal subtraction is critical for gain reconstruction  
- SPE calibration defines the absolute charge scale  
- Fit stability is ensured via constrained initialization  
- The pipeline is modular and reproducible

---

## Author

Developed and maintained by Erasmo Pacini
