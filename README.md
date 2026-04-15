## Usage
### Cloning the Repository 

``` bash
git clone https://github.com/lvc-core/wavelet-analysis.git
```

### Compiling the source code

``` bash
make MODIFICATION NEEDED! (filename should be given via terminal)
```

### Plotting (gnuplot):

``` bash
plot 'data/results_FD.dat' using 1:($2*1000):3 with image
```

# Wavelet analysis 
The wavelet analysis is used in order to get a time resolved frequency spectrum of a signal. Note that the Fourier-Transformation also yields the frequency spectrum of a signal, however this spectrum is not time resolved.    
The idea is that you create a `wavelet` (which is basically just a wavepacket)  and calculate the convolution with the signal at all times $t$ for varying frequencies of the wavelet.

The wavelet analysis has applications in many fields whereas the most common ones are probably signal-processing and image manipulation.  
## Theory
Given any  
* `signal` $f(t)$  
* specific `filter` $\Psi(n)$ with scale $s$ centered at $t_0$  

The convolution of those signals is given by:

$$
C(t_0, s) = \frac{1}{\sqrt{s}} \int_{-\infty}^{\infty} \Psi^* \biggl( \frac{t - t_0}{s} \biggr) f(t) ~dt
$$

A typical choice for the wavelet in physics is the so called `Morlet wavelet` which has the following form:

$$
\Psi^{mother}(t) = \frac{1}{\pi^{1/4} \sqrt{d}} \hspace{0.15cm} \mathrm{exp} \biggl(- \frac{t^2}{2d^2} \biggl) \hspace{0.15cm} e^{i \omega_0 t}
$$

The `daughter wavelet` then describes a scaling and shift in time such that $\int |\Psi(t)|^2 ~dt = 1$

$$
\Psi^{daughter}_{m, n}(t) = \frac{1}{\sqrt{s_m}} \hspace{0.15cm} \Psi^{mother} \biggl( \frac{t - t_n}{s_m} \biggr)
$$

Now it is obvious that the integrand of the convolution is given by the product of the signal and the complex conjugate of a daughter wavelet.

However this calculation is pretty heavy in terms of resources since an integral has to be calculated for each scale.  
This integral can be reduced to a product in frequency space after a Fourier Transformation with respect to $t_0$:

$$
\tilde{C}(\omega, s) = \frac{1}{\sqrt{2\pi}} \int C(t_0, s) ~e^{-i \omega t_0} ~dt_0 = \dots = \sqrt{2 \pi s}  ~ \tilde{\Psi}(s \omega) ~ \tilde{f}(\omega)
$$

# Overlap-save method
Datasets in science can get pretty large and hundres of GB or a few TB are really common. So here it is not guaranteed that the data necessarily fits into the RAM when performing calculations. In order to address this problem the so called Overlap-save method is introduced.

When using the Overlap-save method the data is divided into chunks of a suitable size. Then a wavelet analysis is performed on each chunk individually. Afterwards the inidividual chunks can't simply be concatenated since this would create false data because the edge effects of the chunks are not taken into account. Therefore the procedure is the following:

1. A Morlet-wavelet of length `M` with a specific scale is created. The range of the used scales depends on the frequency-range which should be analyzed.  

2. Calculate the chunk size `N` which is typically a power of 2 so the FFT can use a faster algorithm. The used formula to calculate N is `N = L + (M_max - 1)`. The dataset is divided into blocks of size `L` which is also the number of usable data for each chunk. 
So in other words the data is divided into blocks of size L. For each block the previous M-1 datapoints are also used. The L + M-1 datapoints then make up a chunk. 
Note: For each scale s the block length L varies so for different scales there is a different number of usable data points.
You can find a very good visual explanation of the Overlap-save method [here](https://blog.robertelder.org/overlap-add-overlap-save/).  

3. For each chunk a wavelet analysis is performed and only the last L elements of the result are appended to the result.




## Resources

General overview: [wikipedia](https://en.wikipedia.org/wiki/Overlap%E2%80%93save_method)  
Visualization: [interactive website](https://blog.robertelder.org/overlap-add-overlap-save/)  
Theory: [textbook](https://link.springer.com/book/10.1007/978-3-319-61088-7) chapter 8.4  







