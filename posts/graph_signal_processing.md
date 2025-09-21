---
layout: page
title: Graph Signal Processing
permalink: /posts/graph_signal_processing/
hide: true
---

## Overview

Graphs provide a powerful and intuitive framework for representing relationships in data. Traditionally, signal processing methods 
such as Fourier transform have been used for time series data. But noisy signals are not limited to time series data. 
In this post, I will introduce signal processing on graphs, which extends classical signal processing to data defined on graphs.

* I will first define the Laplacian on graphs. Graph Laplacian is often introduced in a very unintuitive manner via
formulas. But by first constructing graph analogues of gradient and divergence, there's actually a very intuitive derivation. 
* I will then briefly introduce the classical Fourier transform and show how the graph Laplacian provides a similar
framework for signal processing on data defined on graphs. 
* Finally, I will present two potential applications: **images** and **spatial transcriptomics**. 

## Calculus on Graphs

Signal processing methods such as the Fourier transform are defined on the real line, relying on tools from analysis. 
To develop equivalent methods for data defined on graph structures, we first need to establish analogous concepts on graphs. 

### Gradient

In multivariate calculus, each component of the gradient represents the change in function value towards the direction. 
In graphs with associated signals on each node, the equivalent definition can be the change in signal value across edges. 

Formally, we consider an undirected, weighed graph  \\(G = (V, E)\\) where each edge \\(e = (i, j)\\) has some weight \\(W_{ij} > 0\\).

Let \\(f:V\rightarrow \mathbb{R}\\) be a signal on the nodes. We can define the derivative with respect to edge \\(e = (i, j)\\) as:
\\[
 \left(\frac{\partial f}{\partial e}\right) := \sqrt{W_{ij}}(f(j) - f(i))
\\]

Intuitively, a large weight between nodes implies a stronger similarity, so the difference between signals is emphasized (just like in algebra, 
slope is computed as change in y over x, so a smaller change in x means higher slope). 

The gradient operator is the collection of edge derivatives over all edges:
\\[
    \nabla f = \left[ \frac{\partial f}{\partial e_1}, \ldots, \frac{\partial f}{\partial e_m} \right]
\\]
where \\(m\\) is the number of edges in the graph. 

### Divergence

Divergence measures the net outgoing flow of a function from a point. For a graph, that might mean the sum of
the outgoing flow from a node to its neighboring edges. 

With this intuition in mind, for node \\(i\\) in a graph, we define the divergence as:

\\[
    \mathrm{div}\_i(\nabla f) = \sum\_{j \in \mathcal{N}\_i} \sqrt{W\_{ij}} (\nabla_{(i, j)}f)
\\]



where \\(\mathcal{N}\_i\\) is the set of nodes which have an edge to node \\(i\\), and \\(\nabla\_{(i, j)}\\) is the derivative along edge \\(e  = (i, j)\\). 

Put simply, divergence is defined as the sum of the derivatives along all edges connected to the node, scaled by the squared weights of each edge. This aligns with
our intuition that divergence is the net "outgoing flow" from a point.

### Laplacian

Typically, the Laplacian operator is defined as the divergence of the gradient:
\\[
    \Delta := \mathrm{div}\circ\nabla
\\]
So for a graph with signal \\(f\\), we obtain the graph Laplacian as:

$$
\begin{align}
    (\Delta f)(i) &= \sum_{j \in \mathcal{N}_i} \sqrt{W_{ij}} (\nabla_{(i, j)}f) \\
                    &= \sum_{j \in \mathcal{N}_i} W_{ij}(f(i) - f(j))
\end{align}
$$

Intuitively, the graph Laplacian gives us a sense of how the signals at each node differ from the neighboring nodes (i.e. smoothness of signal
at the node). 

### Deriving the Matrix Form of the Laplacian

Our definition of Laplacian, presented above, allows for a nice matrix representation. 

$$
\begin{align}
(\Delta f)(i) &= \sum_{j \in \mathcal{N}_i} W_{ij}(f(i) - f(j)) \\
&= f(i) \sum_{j \in \mathcal{N}_i} W_{ij} - \sum_{j \in \mathcal{N}_i} W_{ij} f(j) \\
&= D_{ii} f(i) - \sum_{j \in \mathcal{N}_i} W_{ij} f(j)
\end{align}
$$

Here, $D$ is the **degree matrix**, a diagonal matrix with entries $D_{ii} = \sum\_{j \in \mathcal{N}\_i} W\_{ij}$, and $W$ is the weighted adjacency matrix. Therefore,

$$
(\Delta f)(i) = (Df)(i) - (Wf)(i)
$$

So in matrix form:

$$
\Delta = D - W
$$

### Positive Semi-Definiteness of the Laplacian
As we will see later, a matrix being positive semi-definite induces a number of convenient properties, such as the eigenvalues being non-negative and the existence of [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) (think of it like a square root of a matrix -- we can only take square roots of positive numbers, unless we consider complex numbers). We now prove that the Laplacian matrix $ L = D - W $ is positive semi-definite. That is, :

$$
\forall f \in \mathbb{R}^n, \quad f^\top L f \geq 0
$$

Clearly, $ L $ is symmetric because both $ D $ and $ W $ are symmetric (the graph is undirected). Now compute:

$$
\begin{align}
f^\top L f &= f^\top (D - W) f \\
&= f^\top D f - f^\top W f
\end{align}
$$

Expanding each term:

$$
\begin{align*}
f^\top D f &= \sum_i D_{ii} f(i)^2 = \sum_{i,j} W_{ij} f(i)^2 \\
f^\top W f &= \sum_{i,j} f(i) W_{ij} f(j)
\end{align*}
$$

Taking the difference:

$$
\begin{align*}
f^\top L f &= \sum_{i,j} W_{ij} f(i)^2 - \sum_{i,j} W_{ij} f(i) f(j) \\
&= \sum_{i,j} W_{ij} \left( f(i)^2 -f(i)f(j) \right)
\end{align*}
$$

Since $W$ is symmetric:

$$
\begin{equation}
    \sum_{i,j} W_{ij} \left( f(i)^2 -f(i)f(j) \right) = \sum_{i,j} W_{ji} \left( f(j)^2 -f(j)f(i) \right)
\end{equation}
$$

Thus: 

$$
\begin{align*}
f^\top L f &= \frac{1}{2} \sum_{i,j} W_{ij} \left( f(i)^2 + f(j)^2 - 2f(i)f(j) \right) \\
&= \frac{1}{2} \sum_{i,j} W_{ij} (f(i) - f(j))^2
\end{align*}
$$

which is clearly always non-negative. Thus, $ L $ is positive semi-definite. 

The equation derived here, $f^TLf = \frac{1}{2} \sum_{i,j} W_{ij} (f(i) - f(j))^2$, is referred to as the discrete p-Dirichlet form of f when p = 2. 

Upon inspection, we notice that it is equal to 0 if and only if the signals are equivalent on all connected nodes. More importantly, the expression is low if the signals of neighboring nodes are similar, and thus encodes the overall smoothness of the graph signals. 

Since $L$ is positive semi definite, there exists a matrix factorization such that $L = KK^T$. 

### Courant-Fisher Theorem

Using the [Courant-Fischer Theorem](https://en.wikipedia.org/wiki/Min-max_theorem), the eigenvalues $ \lambda_0, \lambda_1, \ldots, \lambda_{N-1} $, ordered in **increasing** fashion, of the eigenvalues and the eigenvectors of the graph Laplacian interpreted as follows:

$$
\begin{align}
\lambda_0 &= \min_{f \in \mathbb{R}^N, \|f\|_2 = 1} f^\top L f, \\
\lambda_l &= \min_{\substack{f \in \mathbb{R}^N \\ \|f\|_2 = 1 \\ f \perp \mathrm{span}\{u_0, \ldots, u_{l-1}\}}} f^\top L f, \quad l = 1, 2, \ldots, N-1
\end{align}
$$

where the respective minimums are achieved via $ u_l $, the eigenvector corresponding to eigenvalue $ \lambda_l $.

This theorem, together with the notion that $f^TLf$ represents overall smoothness of the graph signal, gives us an intuitive way to decompose the graph signals. Specifically, eigenvectors corresponding to high eigenvalues signify directions of the signal that are "high frequency". 

### Signal Processing via Fourier Transform. 

Before we dive into signal processing on graphs, I want to first do a quick overview on "classical" signal processing via Fourier transform. 

In classical signal processing, a signal is decomposed into frequencies using the Fourier transform. The basis functions, which represent different frequencies, are complex exponentials:

$$
\begin{equation}
\{e^{2\pi i \xi t}\}_{\xi\in \mathbb{R}}
\end{equation}
$$

These are eigenfunction of differential operators, particularly the 1D Laplacian (second derivative):

$$
\begin{align*}
\Delta e^{2\pi i \xi t} &= \frac{\partial^2}{\partial t^2} e^{2\pi i \xi t} \\
&= -(2\pi \xi)^2 e^{2\pi i \xi t}
\end{align*}
$$

The Fourier transform for each frequency, $\xi$, is defined as:

$$
\begin{align*}
\hat{f}(\xi) &= \langle f, e^{2\pi i \xi t}\rangle\\
&= \int_{\mathbb{R}} f(t) e^{-2\pi i \xi t} dt
\end{align*}
$$

signifying the strength of the signal with the frequency $\xi$. Given the decomposed signal, we can reconstruct the original signal via the inverse Fourier transform :

$$
\begin{align*}
f(t) &= \int_{\mathbb{R}}\langle f, e^{2\pi i \xi t}\rangle e^{2\pi i \xi t}d\xi\\
&= \int_{\mathbb{R}} \hat{f}(\xi) e^{2\pi i \xi t} d\xi
\end{align*}
$$

Most "real" signals such as speech, music, and images, change **gradually** over time or space. This smoothness means their Fourier transform is concentrated at low frequencies. However, signals originating from noise occur irrespective of the underlying "true" signal, and therefore tend to exhibit high frequency. Thus, for the purpose of filtering signals, suppressing signals which exhibit high frequency is of interest. 

In classical filtering, we multiply the signal in the frequency domain by a transfer function:

$$
\begin{equation}
\hat{f}_{\text{out}}(\xi) = \hat{f}_{\text{in}}(\xi) \cdot h(\xi)
\end{equation}
$$

Here,

* $ h(\xi) $ is the transfer function (the "filter")
* $ \hat{f}_{\text{in}} $ is the original signal
* $ \hat{f}_{\text{out}} $ is the processed signal

If $h(\xi)$ dampens high frequencies (e.g. $ h(\xi) \rightarrow 0 $ as $ \|\xi\| \rightarrow \infty$), then high-frequency components are suppressed in the output. This is the basis for designing filters such as low-pass filters in signal processing. To apply a similar filtering process to graphs, we would like to construct a graph analogue of this framework.

### Fourier Transform on Graphs

Now, we want to apply the idea of Fourier transform to graphs. 

Previously, we introduced the graph Laplacian as the Laplacian operator on graphs. In a continuous, Eulicidan space, we observed that the basis for Fourier transform were the eigenfunctions of the Laplacian operator, with their eigenvalues signifying the respective frequency. We aim to construct a similar framework for graphs using the graph Laplacian.  

In analogy with the classical Fourier transform, we interpret the eigenvalues and eigenvectors of \( L \) as generalizations of frequencies and Fourier modes.

Let $ L $ have the eigendecomposition:

$$
L = U \Lambda U^\top,
$$

where $ U = [u_1, \ldots, u_n] $ is an orthonormal matrix of eigenvectors and $ \Lambda = \mathrm{diag}(\lambda_1, \ldots, \lambda_n) $ is the diagonal matrix of eigenvalues. 

We can interpret each eigenvalue $ \lambda_i $ as a “graph frequency,” where larger $ \lambda $ corresponds to higher-frequency components—i.e., more rapid variation over the graph structure.

Any signal $ f \in \mathbb{R}^n $ can be expanded in this eigenbasis:

$$
f(i) = \sum_{k=1}^{n} \hat{f}(k) \, u_k(i),
$$

where

$$
\hat{f}(k) = \langle f, u_k \rangle = \sum_{i=1}^{n} f(i) \, u_k(i)
$$

is the graph Fourier transform corresponding to the "frequency" $ \lambda_k $.

Similar to the classical signal processing as described in the previous section, to filter a graph signal, we apply a transfer function $ h(\lambda) $ which decays as $\lambda\rightarrow \infty$ to each component:

$$
\hat{f}_{\text{out}}(k) = h(\lambda_k) \, \hat{f}(k).
$$

The filtered signal is then reconstructed as:

$$
f_{\text{out}}(i) = \sum_{k=1}^{n} \hat{f}(k) \, h(\lambda_k) \, u_k(i).
$$

For example, such transfer function $h(\lambda)$—which suppresses high-frequency components—can be defined via:

$$
h(\lambda) = \frac{1}{1 + \tau \lambda},
$$

where $ \tau > 0 $ is a tunable parameter that controls the degree of smoothing.

This framework allows us to process signals on arbitrary graph domains in a way that parallels classical Fourier analysis. 

### Application to Image Denoising

Graphs provide a natural framework for representing many types of datasets, and have been widely used in particular for images. In this work, I apply techniques from graph signal processing to denoise images by leveraging information from neighboring pixels.

As a proof of concept, I use the classic “Barbara” image.

```python
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage import data, img_as_float
from skimage.util import random_noise
from skimage.color import rgb2gray
from scipy.ndimage import gaussian_filter
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
from pygsp import graphs, filters
from scipy.sparse.linalg import spsolve
from scipy.sparse import identity

# load image
data = imread('data/barbara.jpg', as_gray=True)
image = img_as_float(data)
plt.imshow(image, cmap='gray')
```

![Barbara Original](/public/images/gsp/barbara.png)

To simulate the denoising procedure, I add random Gaussian noise to each pixel value and try to recover the original image:

```python
# noisy image
noisy_image = random_noise(image, mode = 'gaussian', var = 0.005)
plt.imshow(noisy_image, cmap='gray')
```
![Barbara Noisy](/public/images/gsp/barbara_noisy.png)

One of the most common methods used to denoise images is the Gaussian filter, where each pixel is essentially a weighed average of neighboring pixels using the 2D Gaussian kernel. 

```python
# gaussian filter
filtered_image = gaussian_filter(noisy_image, sigma=3)
plt.imshow(filtered_image, cmap='gray')
```
![Barbara Gaussian](/public/images/gsp/barbara_gaussian.png)

While Gaussian filter does 'smooth out' the noise, it tends to blur the image excessively (this can be tuned via adjusting sigma). Moreover, there still exists
noticeable texture from the Gaussian noise. 

Now, we implement the graph signal processing approach described above. In the context of images, each pixel represents a node, and neighboring nodes are connected via an edge. 
We define the relevant functions below:
```python
# construct graph
def img_to_graph_default(img, theta = 0.1):
    rows, cols = img.shape
    graph = lil_matrix((rows * cols, rows * cols))

    # define neighbors
    neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1),  
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]  
    for r in range(rows):
        for c in range(cols):
            idx = r * cols + c
            for dr, dc in neighbors:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    nidx = nr * cols + nc
                    graph[idx, nidx] = 1

    return graph.tocsr()

def tikhonov_denoising_default(image, alpha):
    
    y =  image.flatten()
    rows, cols = image.shape
    W = img_to_graph_default(image)
    G = graphs.Graph(W)  
    G.compute_laplacian('combinatorial') 
    L = G.L 
    
    I = identity(L.shape[0], format='csr')
    A = I + alpha * L
    x_hat = spsolve(A, y)
    
    return x_hat.reshape(rows, cols)
```
Let us now run the denoising process on the Barbara image
```python
hoge = tikhonov_denoising_default(image, 10)
plt.imshow(hoge,cmap='gray')
```
![Barabara Tikhonov Default](/public/images/gsp/barabara_tikhonov_default.png)

Indeed, the image appears smoother after denoising. However, it is clear that we have **over-denoised**—the result looks overly blurred and out of focus. This blurring is a common obstacle when working on signal imputation. Most imputation methods rely on "borrowing" information from neighboring elements, typically through some form of weighted averaging. While this reduces noise, it also inevitably smooths out fine details, leading to a loss of granularity and, consequently, important information. 

So how do we work around this?

One solution I found to be particularly interesting was to **rely only on neighboring signals that are similar in value to the original signal, rather than averaging across all nearby elements.**

Consider image denoising: the main reason images become blurry is the loss of edge features, which are sudden changes in pixel values between neighbors. This makes sense—our graph signal processing formulation penalizes sharp changes between neighboring pixels, so traditional algorithms tend to smooth edges away.

But if we assume that noise is roughly Gaussian, then large jumps between pixel values are unlikely to come from noise. Instead, noise usually appears as smaller, moderate fluctuations. That means that large fluctuations in signal are often real signal, not artifacts.

So instead of averaging across all nearby pixels, we should average only among pixels that are similar in value. This way, noise is reduced while edges are preserved.

<img src="/public/images/gsp/gsp-anchor.jpeg" alt="GSP Anchor" width="300"/>

In the figure above, suppose we want to impute the value of the center pixel. If we average the three pixels on the right, the edge gets blurred. But if we average the remaining neighbors—whose values are closer to the center pixel—we preserve the edge while still reducing variance via weighted averaging of data points.

We can implement this solution as follows, where we weigh each edge by the pixel similarity:

```python

# compute edge weight
def compute_weight(x, y, theta=0.1):
    difference = np.abs(x-y)
    if difference >= 0.05:
        return 0
    else:
        return np.exp(-(difference**2)/(2*theta**2))

# construct graph
def img_to_graph_custom(img, theta = 0.1):
    rows, cols = img.shape
    graph = lil_matrix((rows * cols, rows * cols))

    # define neighbors
    neighbors = [(-1, 0), (1, 0), (0, -1), (0, 1),  
                 (-1, -1), (-1, 1), (1, -1), (1, 1)]  
    for r in range(rows):
        for c in range(cols):
            idx = r * cols + c
            for dr, dc in neighbors:
                nr, nc = r + dr, c + dc
                if 0 <= nr < rows and 0 <= nc < cols:
                    nidx = nr * cols + nc
                    weight = compute_weight(img[r, c],img[nr, nc])
                    graph[idx, nidx] = weight

    return graph.tocsr()


def tikhonov_denoising(image, alpha):
    
    y =  image.flatten()
    rows, cols = image.shape
    W = img_to_graph_custom(image)
    G = graphs.Graph(W)  
    G.compute_laplacian('combinatorial') 
    L = G.L 
    
    I = identity(L.shape[0], format='csr')
    A = I + alpha * L
    x_hat = spsolve(A, y)
    
    return x_hat.reshape(rows, cols)
```

Applying this to our noisy image:

```python
hoge = tikhonov_denoising(image, 10)
plt.imshow(hoge,cmap='gray')
```
![barbara tikhonov](/public/images/gsp/barbara_tikhonov_alt.png)

Now we see that the image is significantly smoother, yet it still preserves important features. For example, the stripes in the tablecloth, scarf, and pants remain clearly visible.

### Application to Spatial Transcriptomics

Another area in which there have been extensive work on signal imputation is single-cell and spatial omics data (e.g. [MAGIC](https://www.cell.com/cell/fulltext/S0092-8674(18)30724-4), [scImpute](https://www.nature.com/articles/s41467-018-03405-7)). Biological data are often noisy. In particular,
single-cell and spatial transcriptomics data are notorious for their noise, most notably from dropout. I've worked extensively with these datasets, so I was curious to see if the graph signal processing approach
could be effective. 

Here, we use the pre-processed 10x Genomics Visium H&E dataset of adult mouse brain. The raw data can be downloaded [here](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain?)

```python
import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import squidpy as sq
import scipy
from pygsp import graphs, filters
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")
```
```python
# load the pre-processed dataset
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
```
Let us first look at the raw signal for several marker genes:
```python
sq.pl.spatial_scatter(adata, color=["Reln", "Cux1", "Necab1", "Foxp2"])
```
<img src="/public/images/gsp/ST_raw.png" alt="ST Raw" width="1500"/>

While we do observe certain spatial expression patterns, the signals are very noisy. This is because spatial transcriptomics experiments are inherently sparse and noisy due to technical and biological factors: only a fraction of transcripts in each cell are captured, and sequencing depth is limited.

To denoise this data, we could employ a similar strategy as we did with the Barbara image.
First, obtain spatial coordinates of each data point, and look for **spatial neighbors** e.g. pixels that are close in terms of spatial location:

```python
# get PC's
spatial_location = adata.obsm['spatial']
# compute distances between spatial locations
dist = scipy.spatial.distance_matrix(spatial_location, spatial_location)
# sort by distance and get coordinates
nbrs_sorted = np.argsort(dist,axis=1)
```
Then, for each pixel, we compute **expression distance** from every other pixel e.g. how similar their expression patterns are across the transcriptome. 
To avoid the curse of dimensionality, we compute the Euclidean distance on PCA space.
```python
# compute distance in GEX space
gex = adata.obsm['X_pca']
dist_gex = scipy.spatial.distance_matrix(gex, gex)
```
Using this information, we can construct the graph as follows and perform denoising:
```python
# construct graph
cell_num = adata.n_obs
graph = scipy.sparse.lil_matrix((cell_num, cell_num))

# compute sigmas
sigmas = np.zeros(cell_num)
for idx in range(cell_num):
    spatial_nbrs = nbrs_sorted[idx, 1:7]
    nbrs_distances = dist_gex[idx, spatial_nbrs]
    sigmas[idx] = np.sort(nbrs_distances)[-1]

# construct graph
for idx in range(cell_num):
    spatial_nbrs = nbrs_sorted[idx, 1:7]
    for nbr in spatial_nbrs:
        dist = dist_gex[idx, nbr]
        adaptive_sigma = sigmas[idx]*sigmas[nbr]
        if adaptive_sigma == 0:
            d = 0
        else: 
            d = np.exp(-1*dist**2/adaptive_sigma)
        graph[idx, nbr] = d
        graph[nbr, idx] = d

G = graphs.Graph(graph)  
G.compute_laplacian('combinatorial') 
L = G.L 
I = scipy.sparse.identity(L.shape[0], format='csr')
A = I + 10 * L

denoised = scipy.sparse.linalg.spsolve(A, adata.X)
adata.layers['denoised'] = denoised
```
Here, we employ a method for adaptive tuning of the hyperparameter sigma, presented in the paper [Self-Tuning Spectral Clustering](https://papers.nips.cc/paper_files/paper/2004/file/40173ea48d9567f1f393b20c855bb40b-Paper.pdf), to handle heterogeneous variability across cell types. 

We can now visualize the denoised expression results:
```python
sq.pl.spatial_scatter(adata, color=["Reln", "Cux1", "Necab1", "Foxp2"], ncols=2)
```
<img src="/public/images/gsp/ST_denoised.png" alt="ST denoised" width="1500"/>

### Discussion

We see that our graph signal processing (GSP) algorithm does a nice job of denoising spatial transcriptomics data. The spatial expression patterns of each gene look much clearer and smoother now.

I think this makes for an interesting proof-of-concept of how GSP could be applied in computational biology. That said, we’ll definitely need a stronger benchmark to really test whether this approach improves spatial transcriptomics analysis—and whether it does so meaningfully better than existing imputation methods.

It’s worth noting that several studies have found that scRNA-seq imputation often doesn’t significantly improve downstream analysis and can even lose important biological signals (for example, this [one](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x)
).

On top of that, the current GSP algorithm we applied to spatial data uses some extra tricks—like edge weighting and self-tuning spectral clustering—to boost performance. In theory, those additions should help the model capture the data structure better. But to be sure, it’d be really useful to run an ablation study to test whether those components actually make a difference.

Feel free to reach out if you have any questions or ideas—I’d love to hear your thoughts!