This software contains the implementation of the multi-view kernel k-means (MVKKM) and multi-view spectral clustering (MVSpec)
algorithms described in G.Tzortzis and A.Likas, "Kernel-based Weighted Multi-view Clustering", ICDM 2012.

The main routines for these two algorithms can be found in the MVKKM.m and MVSpec.m files.
Type "help MVKKM" or "help MVSpec" within MatLab for usage details.

Moreover, a demo file (Demo.m) is included that shows how to call the above functions.
The synthetic dataset of our paper is used for the demonstartion.
It contains two views and three clusters (synthetic_data.mat).
Also, a precomputed rbf kernel, with sigma=0.2, is provided for each view (synthetic_data_kernel.mat).

The code for the global kernel k-means algorithm, that was used to initialize MVKKM
in our experiments, is also provided (Weighted_Global_Kernel_K_Means.m).
For more details on this method please refer to 
G.Tzortzis and A.Likas, "The Global Kernel k-Means Algorithm for Clustering in Feature Space", IEEE TNN, 2009.

Finally, the kernel k-means and global kernel k-means routines are taken directly from our previous software,
which can also be found at: http://www.iit.demokritos.gr/~gtzortzi.

For any questions regarding the implementation do not hesitate to contact me at: gtzortzi@iit.demokritos.gr

G. Tzortzis