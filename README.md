# 🌈 Bayesian Fully Connected Tensor Network for Hyperspectral and Multispectral Image Fusion

**Authors:** Linsong Shan, Zecan Yang, Laurence T. Yang, Changlong Li, Honglu Zhao, Xin Nie  
📄 **Paper:** [Bayesian Fully Connected Tensor Network for Hyperspectral and Multispectral Image Fusion](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=11540397)  


---

## 🧭 Introduction
<p align="center">
  <img src="https://github.com/user-attachments/assets/a080f005-d221-461b-ac33-fb70acf80dec"
       alt="image"
       width="50%">
</p>

**Bayesian Fully Connected Tensor Network (BFCTN)** is a novel Bayesian framework designed for **hyperspectral–multispectral image fusion**.  
Unlike traditional tensor-based or deep learning fusion models, BFCTN introduces a **fully connected tensor representation** that jointly models spectral, spatial, and latent dependencies under a **Bayesian inference paradigm**.

### 🔍 Key Highlights
- **Bayesian Tensor Modeling**: Incorporates uncertainty modeling and adaptive regularization through Bayesian inference.  
- **Fully Connected Structure**: Establishes inter-layer correlations across all tensor modes, enabling global information propagation.   
- **Efficient Optimization**: Employs closed-form updates for posterior estimation with low computational overhead.


## 🗂 Folder Structure

```plaintext
BFCTN/
├── Model/              # Model implementation
├── Demo.m              # A simple demo to test the method
├── Function/           # Utility functions
├── results/            # Fusion results
└── README.md           # This file
```

## 🚀 Getting Started

To reproduce the basic experiment and fusion results, run:

```matlab
Demo
```

## 📦 Requirements

MATLAB R2022a or later

## ✅ Key Features

- Implementation of the Bayesian Fully Connected Tensor Network (BFCTN) for hyperspectral and multispectral image fusion.
- Support for various spatial degradation models and noise configurations.
- Modular, easy-to-extend code for testing new priors, tensor decompositions, or fusion datasets.
- Reproducible and interpretable Bayesian estimation steps.


## 📊 Example Results
<details>
<summary>🧪 Click to Expand Example Results</summary>

<img width="1238" height="691" alt="image" src="https://github.com/user-attachments/assets/e45ac3fe-b5f4-4c66-8971-7a2f0dc49d95" />

<img width="1210" height="683" alt="image" src="https://github.com/user-attachments/assets/3ceacc79-875c-400c-b7d5-d45c0cab68b7" />


*(Refer to the paper for full benchmark results and analysis.)*

</details>

## 🔬 Citation

If you find this work helpful, please consider citing our paper:
```
@article{shan2026bayesian,
  title={Bayesian Fully-Connected Tensor Network for Hyperspectral-Multispectral Image Fusion},
  author={Shan, Linsong and Yang, Zecan and Yang, Laurence T and Li, Changlong and Zhao, Honglu and Nie, Xin},
  journal={IEEE Transactions on Image Processing},
  year={2026},
  publisher={IEEE}
}
```


