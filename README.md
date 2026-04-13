# tesi-pmt-gain_characterization



This project implements a full analysis pipeline for the characterization of Photomultiplier Tubes (PMTs) used in high-energy physics experiments.

The goal is to extract the PMT gain as a function of applied high voltage (HV) through a complete calibration workflow including pedestal subtraction, single-photoelectron (SPE) calibration, and gain curve reconstruction.

The analysis is based on ROOT and follows standard techniques used in detector calibration in astroparticle physics.

---

## Data

The dataset consists of experimental PMT acquisition runs under different operating conditions:

- **Pedestal runs** (no light conditions)  
  Used to estimate the electronic baseline (ADC offset) for each channel.

- **Single-photoelectron (SPE) runs**  
  Low-intensity light data used to calibrate the charge response in terms of photoelectrons.

- **Gain curve runs (high illumination)**  
  High-intensity data acquired at different high voltages (900–1500 V) to reconstruct the PMT gain curve.

All raw data are stored in CSV format and processed using ROOT-based histograms.

---

## Requirements

The following software and Python packages are required:

- ROOT (v6 or higher)
- numpy
- pandas
- pathlib
- IPython (for Jupyter/Colab visualization)

---

## Analysis Pipeline

The analysis is structured into three main stages:

### 1. Pedestal Estimation

- Extraction of electronic baseline per channel  
- Gaussian fit of ADC distributions under dark conditions  
- Determination of mean pedestal value (μ_ped)

---

### 2. Single-Photoelectron (SPE) Calibration

- Pedestal subtraction from low-light data  
- Modeling of charge spectrum using a double-Gaussian function  
- Extraction of SPE mean charge per channel  
- Conversion from ADC units to physical charge (pC)

---

### 3. Gain Curve Reconstruction

- Analysis of charge distributions as a function of high voltage  
- Gaussian fit of high-intensity charge spectra  
- Conversion to number of photoelectrons using SPE calibration  
- Computation of PMT gain:

$$
G(V) = \frac{Q(V)}{N_{pe} \cdot e}
$$

- Empirical fit of gain curve using a power-law model:

$$
G(V) = A \cdot V^ b
$$

- Extraction of characteristic HV points for each channel

---

## Physical Quantities Extracted

The pipeline allows estimation of:

- Electronic pedestal (ADC baseline)
- Single-photoelectron charge response
- Number of photoelectrons per channel
- PMT gain as a function of HV
- Power-law scaling exponent of gain curves
- Characteristic operating voltage per channel

---

## Repository Structure
