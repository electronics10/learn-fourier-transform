# learn-fourier-transform
Learn Fourier transform by doing!

## Continuous Fourier transform

By definition, $$\hat{x}(f)=\int_{\mathbb{R}}x(t)e^{-i2\pi ft}dt.$$ For example, the Fourier transform of a sine function $\sin(2\pi f_0 t)$ is $$\int_{\mathbb{R}}\frac{e^{i2\pi f_0t}-e^{-i2\pi f_0t}}{2i}e^{-i2\pi ft}dt=\frac{1}{2i}\int_{\mathbb{R}}[e^{-i2\pi (f-f_0)t}-e^{-i2\pi (f+f_0)t}]dt=\frac{\delta(f-f_0)-\delta(f+f_0)}{2i}.$$ (The integration of exponential isn't an ordinary convergent integral, it's simply a distributional identity to the delta functional; accept it for now.)

  

Note that the field we operate in is the complex number system. Since $t$ and $f$ are mostly real, it is useful to picture the functions in the time domain to be a combination of spirals ($e^{i2\pi ft}$); and picture the functions in the frequency domain also as an 3D image where each spike at a specific frequency is the value/density of the spiral of that frequency $e^{i2\pi ft}$ while the complex value indicate the phase, e.g. spike at frequency $f_0$ with value $2i$ means that it produces a spiral $2ie^{i2\pi ft} = 2e^{i(2\pi ft+\frac{\pi}{2})}$. ($e^{i2\pi ft}$ provide $\cos$ in the real part by default.)

  

Physically, the result $\dfrac{\delta(f-f_0)-\delta(f+f_0)}{2i}$ says that a pure tone of frequency $f_0$ contains energy at exactly two frequencies, $+f_0$ and $-f_0$, and nothing else. Each delta is an infinitely sharp, infinitely tall spike of zero width — the spectral signature of a signal that oscillates forever at a single rate.

  

Three features deserve unpacking, because each will reappear (in modified form) the moment we move from this idealized integral to a real computation.

  

**The two spikes, and why frequency comes in pairs.** A real-valued signal like $\sin(2\pi f_0 t)$ always has a spectrum with this mirror symmetry: a peak at $+f_0$ is accompanied by a partner at $-f_0$. The reason traces back to Euler's formula. To build a *real* oscillation out of the complex basis $\{e^{i2\pi ft}\}$, you need two complex exponentials spinning in opposite directions — one at $+f_0$, one at $-f_0$ — so that their imaginary parts cancel and only the real sine survives. Negative frequency is not a physical "frequency below zero"; it is the bookkeeping device that lets a complex basis represent real signals. (Hypothesis under which this holds: $x(t)$ is real-valued. For a genuinely complex signal, the symmetry breaks and a single spike can stand alone.)

  

**The factor $\tfrac{1}{2i}$, and why the spectrum is imaginary.** The transform here is purely imaginary and antisymmetric: the coefficient at $+f_0$ is $\tfrac{1}{2i} = -\tfrac{i}{2}$, and at $-f_0$ it is $+\tfrac{i}{2}$. This encodes phase. A sine is a cosine shifted by a quarter period, and that quarter-period shift shows up precisely as the factor of $i$. Had we transformed $\cos(2\pi f_0 t)$ instead, the spikes would be real and equal, $\tfrac{1}{2}(\delta(f-f_0)+\delta(f+f_0))$. So the *location* of the spikes tells you which frequencies are present; the *complex value* at each spike tells you the amplitude and phase.

  

**The delta function, and the trap we are about to fall into.** The clean delta spike is a consequence of integrating over *all* of $\mathbb{R}$ — the signal must oscillate from $t=-\infty$ to $t=+\infty$. This is exactly the assumption a real measurement cannot satisfy. We record for a finite time, at discrete instants, on a finite machine. Each of those three compromises deforms the delta in a specific, predictable way, and understanding those deformations *is* practical Fourier analysis. The rest of this article is the story of how the two ideal spikes turn into what your computer actually plots.

  

## Discrete Fourier transform

The bridge from the ideal integral to the DFT has three planks, laid in a specific order. Each replaces one impossible assumption with a real constraint, and each deforms the spectrum in a named, predictable way.

  

### Plank 1: Finite duration (truncation → leakage)

  

We can't integrate over all of $\mathbb{R}$. We observe $x(t)$ only on $[-T, T)$. Formally, this is multiplying the eternal signal by a **rectangular window**: $$w(t) = \begin{cases} 1 & |t| < T \ 0 & \text{otherwise.}\end{cases}$$ What we actually transform is $x(t)w(t)$, not $x(t)$.

  

According to the comvolution theorem, multiplication in time becomes convolution in frequency: $$\widehat{x \cdot w}(f) = \hat{x}(f) * \hat{w}(f).$$The rectangle's transform is a sinc: $\hat{w}(f) = 2T\,\dfrac{\sin(2\pi T f)}{2\pi T f} = 2T\,\text{sinc}(2Tf)$. Convolving a delta with a function just _places that function at the delta's location_. So each ideal spike becomes a sinc lobe centered at $\pm f_0$: $$\frac{\delta(f-f_0)}{2i} \;\longrightarrow\; \frac{2T}{2i}\,\text{sinc}\big(2T(f-f_0)\big).$$

  

A main lobe of width $\sim 1/T$ at each tone, flanked by ringing **side lobes**. This smearing is spectral leakage: energy that "should" sit at exactly $f_0$ leaks into neighboring frequencies. The narrower the window (small $T$), the wider the lobe — this is the uncertainty principle made visible.

  

### Plank 2: Discrete sampling (sampling → periodic spectrum + aliasing)

  

We don't have a function, we have samples taken every $\Delta t = 1/f_s$ seconds. Sampling is multiplying by a **Dirac comb** (an infinite train of deltas spaced $\Delta t$ apart).

  

We can prove that a comb's transform is another comb, spaced $f_s$ apart in frequency. Multiplying by a comb in time → convolving with a comb in frequency → the spectrum gets copied and repeated at every multiple of $f_s$. So the discrete-time spectrum is periodic with period $f_s$.

  

The hazard: **aliasing.** If those repeated copies overlap — which happens when the signal contains frequencies above $f_s/2$ — they add together and become indistinguishable. A high tone masquerades as a low one. The threshold $f_s/2$ is the Nyquist frequency. Nyquist–Shannon sampling theorem states that if f $x(t)$ contains no frequencies at or above $f_s/2$, then the samples determine $x(t)$ exactly.

  

### Plank 3: Finitely many samples (the DFT itself)

  

Combine Planks 1 and 2: finite and discrete. You have exactly $N$ samples. The only transform you can actually compute is the Discrete Fourier Transform:

  

**Definition** $$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-2\pi i k n/N}\, \qquad k = 0, 1, \dots, N-1.$$This is a finite sum — fully computable. `np.fft.fft` evaluates it (the FFT is just a fast algorithm for the same numbers, status: algorithm not new math).

  

**Two facts you must internalize to read the output.**

  

_Which frequencies $k$ corresponds to._ The DFT only reports values at the discrete grid $$f_k = \frac{k}{N}f_s, \qquad \text{spacing } \Delta f = \frac{f_s}{N} = \frac{1}{N\Delta t}.$$ $\Delta f$ is your **frequency resolution** — and notice it equals $1/(\text{total record length})$, tying directly back to Plank 1's lobe width. Longer recording → finer bins. This is the same uncertainty tradeoff.

  

_The implicit periodicity assumption._ The DFT treats your $N$ samples as **one period of an infinitely repeating signal**. It silently wraps the end back to the start. If your window does _not_ contain a whole number of periods, that wrap creates a discontinuity, and the discontinuity is what produces visible leakage in the discrete spectrum. This is the practical face of Plank 1.

  

---

  

## Example: DFT of a sine wave
### The basics

  Let's take the simplest example where the signal $x(t)=\sin(2\pi f_{0}t)$. We need to first setup time, which is inevitably **finite** and **discrete**.
  ```python
  T = 1.4 # sec, monitored time interval = window
  t = np.linspace(0, T, N+1) # N discrete time segment in [0, T]
  ```

We can obtain the sampling rate accordingly (by the way, `linspace` isn't a good bookkeeping representation.)

```python
T = 1.4 # sec, monitored time interval = window
fs = 200 # sampling rate -> Delta t = 1/fs -> N = T/(Delta t) = T*fs
t = np.arange(0, T, 1/fs) # N discrete time segment in [0, T)
N = len(t) # T*fs may not be an integer
```

Now the signal can be represented as 
```python
f0 = 5 # less than fs/2
x = np.sin(2*np.pi*f0*t)
```

where its DFT is 
```python
X = np.fft.fft(x)
```

Note that it is hard to intepret in a spectral sense since we are actually doing $$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-2\pi i k n/N}\, \qquad k = 0, 1, \dots, N-1.$$
And the result is just a list of $N$ values of coefficient indicators of $e^{2\pi ikt}$. Firstly, $f_k=\frac{k}{N}f_s$ so we need to actually specify the frequency spacing.

```python
freqs = np.fft.fftfreq(N, 1/fs)
```

The DFT is computed and we can show it by

```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2,1, figsize=(9,4))
ax[0].plot(t,x)
ax[0].set_title(f"sin(2pi {f0}t), fs={fs}, T={T}, N={N}")
ax[1].plot(freq, np.abs(X))
ax[1].set_title("Absolute value of the spectrum")

fig.tight_layout()
plt.show()
```

Note that $X$ is complex so we take the absolute value (ignore the phase).

### Nuance
**The governing identity (status: derived from the definition).** The DFT only ever evaluates the sampled spiral $$e^{i2\pi f_k t_n} = e^{i2\pi \frac{kf_s}{N}\frac{nT}{N}} = e^{i2\pi \frac{k}{T}\frac{nT}{N}} = e^{i2\pi k n/N}.$$ As a function of the integer $k$, this is **periodic with period $N$**: $\;k\to k+N$ leaves it unchanged, since $e^{i2\pi n}=1$. Therefore $X[k]=X[k+N]$. The index $k$ does not live on the line segment ${0,\dots,N-1}$; it lives on a **circle**, the finite cyclic group $\mathbb{Z}/N\mathbb{Z}$. Everything below is a corollary.

**Which physical frequency does bin $k$ represent?**

Bin $k$ labels the spiral $e^{i2\pi kn/N}$, nominally at frequency $f_k = \tfrac{k}{N}f_s$. But because $k$ is periodic mod $N$, **bin $k$ is identical to bin $k-N$**, which corresponds to $f_k - f_s$. Each bin therefore stands for _two_ candidate frequencies — one near the top of $[0,f_s)$, one negative near zero.

We resolve the ambiguity by **convention, justified by Nyquist**: read every bin as the frequency of **smallest magnitude**. A real signal sampled at $f_s$ can only carry frequencies in $(-f_s/2,\,f_s/2)$; anything outside has aliased into that band, so the small-magnitude reading is the only physically admissible one. Concretely, for the upper half $k>N/2$: $$k \;\longleftrightarrow\; k-N \quad (\text{a negative frequency}).$$This is exactly why `np.fft.fftfreq(N, 1/fs)` returns the second half as **negative** numbers: it has already performed the $k\to k-N$ relabeling. The array runs $0,\Delta f,\dots,(N/2-1)\Delta f$, then jumps to $-\tfrac{N}{2}\Delta f,\dots,-\Delta f$, with $\Delta f = f_s/N$.

```python
X = np.fft.fft(x)
freqs = np.fft.fftfreq(N, 1/fs)
print(freqs)
print(X)
X_sh  = np.fft.fftshift(X)
f_sh = np.fft.fftshift(freqs)   # apply the SAME shift to the axis
print(f_sh)
print(X_sh)
```

### Inspection
**Aliasing**
```python
import numpy as np
import matplotlib.pyplot as plt

fs = 200
T  = 1.0
t  = np.arange(0, T, 1/fs)
N  = len(t)

f0_list = [5, 95, 105, 195]   # 105 and 195 are ABOVE Nyquist = 100

fig, ax = plt.subplots(1, 4, figsize=(16, 3), sharey=True)
for a, f0 in zip(ax, f0_list):
    x = np.sin(2*np.pi*f0*t)
    X     = np.fft.fftshift(np.fft.fft(x))
    freqs = np.fft.fftshift(np.fft.fftfreq(N, 1/fs))
    a.stem(freqs, np.abs(X)/N)
    f_alias = (f0 + fs/2) % fs - fs/2   # fold into (-fs/2, fs/2]
    a.set(title=f"f0={f0} → reads {f_alias:+.0f}", xlabel="f (Hz)",
          xlim=(-fs/2, fs/2))
ax[0].set_ylabel("|X|/N")
plt.tight_layout(); plt.show()
```

| $f_0$ | $f_0 - f_s$ | reported peak                    |
| ----- | ----------- | -------------------------------- |
| 5     | —           | $\pm 5$ (correct)                |
| 95    | —           | $\pm 95$ (correct)               |
| 105   | $-95$       | $\pm 95$ (looks identical to 95) |
| 195   | $-5$        | $\pm 5$ (looks identical to 5)   |

**Leakage**
```python
import numpy as np
import matplotlib.pyplot as plt

fs, f0 = 200, 5

fig, ax = plt.subplots(1, 3, figsize=(13, 3.2), sharey=True)
for a, T in [(ax[0], 1.0), (ax[1], 1.1), (ax[2], 1.05)]:
    t = np.arange(0, T, 1/fs)
    N = len(t)
    x = np.sin(2*np.pi*f0*t)
    X     = np.fft.fftshift(np.fft.fft(x))
    freqs = np.fft.fftshift(np.fft.fftfreq(N, 1/fs))
    a.stem(freqs, np.abs(X)/N)
    a.set(title=f"T={T} → {f0*T:.1f} periods", xlabel="f (Hz)",
          xlim=(0, 15))
ax[0].set_ylabel("|X|/N")
plt.tight_layout(); plt.show()
```

The ratio $f_0/\Delta f = f_0/(f_s/N) = f_0 \cdot N/f_s = f_0 \cdot T$ collapses to the period count — so "peak lands on a bin," "$f_0/\Delta f$ integer," and "window holds whole periods" are three names for one condition. Stress-test the slogan "non-integer periods → leakage": it is exact for the rectangular window, but the _amount_ of side-lobe energy is a property of the window shape, not a universal constant — a Hann window on the same $T=1.1$ signal would show dramatically lower side lobes (at the cost of a wider main lobe). So the right framing is: _leakage is caused by the implicit periodic wrap creating a discontinuity; window shaping controls how badly that discontinuity rings._
