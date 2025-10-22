# 🌈 Bayesian Fully Connected Tensor Network for Hyperspectral and Multispectral Image Fusion

**Authors:** Linsong Shan, Laurence T. Yang, Zecan Yang, Changlong Li, Honglu Zhao, Xin Nie  
📄 **Paper:** [Bayesian Fully Connected Tensor Network for Hyperspectral and Multispectral Image Fusion](https://arxiv.org/abs/2510.18400)  
[![arXiv](https://img.shields.io/badge/arXiv-2510.18400-b31b1b.svg)](https://arxiv.org/abs/2510.18400)

---

## 🧭 Introduction
<img width="1116" height="748" alt="image" src="https://github.com/user-attachments/assets/a080f005-d221-461b-ac33-fb70acf80dec" />


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
@article{shan2025bfctn,
  title={Bayesian Fully Connected Tensor Network for Hyperspectral and Multispectral Image Fusion},
  author={Linsong Shan, Laurence T. Yang, Zecan Yang, Changlong Li, Honglu Zhao, Xin Nie},
  journal={arXiv preprint arXiv:2510.18400},
  year={2025}
}
```


