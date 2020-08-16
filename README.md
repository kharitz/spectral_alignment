# Graph Spectral Alignment
This repository contains the matlab implementation of the graph spectral alignment originally proposed in "FOCUSR: Feature Oriented Correspondence using Spectral Regularization - A Method for Precise Surface Matching", by Lombaert, H. et. al., published in PAMI 2013. We perform this spectral alignment to overcome the limitationso of spectral graph convolution networks in our work "Graph Convolutions on Spectral Embeddings for Cortical Surface Parcellation", published in Medical Image Analysis, January 2019. 

### What does the reopositery do?
main.m
  1. FreeSurfer brain surfaces is read from the "dataset" folder.
  2. Spectral embedding of the brain graph is computed.
  3. Spectral basis of each subject is aligned to a common reference from the dataset.
  4. Computed embeddings, transformation and aligned spectral embeddings are saved in "output" folder.
mat2pyt.py
  1. Convert data from ".mat" files to ".pt" for DeepLearning algorithms later use in PyTorch.

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

#### REFERENCE 
Please cite our papers if you use this code in your own work:

```
@article{lombaert2012focusr,
  title={FOCUSR: feature oriented correspondence using spectral regularization--a method for precise surface matching},
  author={Lombaert, Herve and Grady, Leo and Polimeni, Jonathan R and Cheriet, Farida},
  journal={IEEE transactions on pattern analysis and machine intelligence},
  volume={35},
  number={9},
  pages={2143--2160},
  year={2012},
  publisher={IEEE}
}
```
```
@article{gopinath2019graph,
  title={Graph convolutions on spectral embeddings for cortical surface parcellation},
  author={Gopinath, Karthik and Desrosiers, Christian and Lombaert, Herve},
  journal={Medical image analysis},
  volume={54},
  pages={297--305},
  year={2019},
  publisher={Elsevier}
}
```

### Future
The python implementation of the graph spectral alignment is coming soon!
