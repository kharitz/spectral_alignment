# Graph Spectral Alignment
This repository contains the matlab implementation of the graph spectral alignment originally proposed in "FOCUSR: Feature Oriented Correspondence using Spectral Regularization - A Method for Precise Surface Matching", by Lombaert, H. et. al., published in PAMI 2013. We perform this spectral alignment to overcome the limitationso of spectral graph convolution networks in our work "Graph Convolutions on Spectral Embeddings for Cortical Surface Parcellation", published in Medical Image Analysis, January 2019. 

### What does the reopositery do?
1. FreeSurfer brain surfaces is read from the "dataset" folder.
2. Spectral embedding of the brain graph is computed.
3. Spectral basis of each subject is aligned to a common reference from the dataset.
4. Computed embeddings, transformation and aligned spectral embeddings are saved in "output" folder.
5. Data from ".mat" files are converted to ".pt" for DeepLearning algorithms in PyTorch.

### Where to find the dataset?
- The MindBoggle dataset is available to download [here](https://osf.io/nhtur/).
- The ADNI dataset is available to download [here](http://adni.loni.ucla.edu).

### Package Requirements
- matlab2018 or higher
- python3, pytorch>1.0 

### Usage
- Open matlab in the downloaded folder and RUN the file "main.m" to perform 1-4.
- Run python script to convert ".mat"-->".pt"
```
python3 mat2pyt.py
```

