"""
LEARNING THE FOURIER TRANSFORM BY DOING
=======================================

Run each section in order. Each one isolates ONE idea.
Every plot is saved as a PNG so you can look at them side by side.

Core object: the Discrete Fourier Transform (DFT)
    X[k] = sum_{n=0}^{N-1} x[n] * exp(-2*pi*i * k*n / N)
This is a CHANGE OF BASIS: it rewrites your N samples in the basis of
N pure complex exponentials. X[k] is the coordinate along the k-th one.
numpy.fft.fft computes exactly this; "FFT" is just the fast algorithm.

Mental model to keep:
    time/space domain  <-- FFT/IFFT -->  frequency domain
    "what the signal is"                 "which oscillations it's made of"
"""

import numpy as np
import matplotlib.pyplot as plt

# A tidy helper so every section produces the same kind of two-panel figure.
def show_time_and_spectrum(t, x, fs, title, fname, real_signal=True):
    """Plot a signal and its magnitude spectrum.

    t  : time samples (seconds)
    x  : signal samples
    fs : sampling rate (samples per second). Sets the frequency axis units.
    """
    N = len(x)
    X = np.fft.fft(x)                 # the DFT, defined above
    freqs = np.fft.fftfreq(N, d=1/fs)  # maps bin index k -> physical frequency (Hz)

    # For a REAL signal the spectrum is conjugate-symmetric, so the negative
    # frequencies are redundant. We show only the non-negative half.
    if real_signal:
        keep = freqs >= 0
        f_plot = freqs[keep]
        mag = np.abs(X[keep]) * 2 / N   # *2/N -> recover real amplitude per bin
    else:
        order = np.argsort(freqs)
        f_plot = freqs[order]
        mag = np.abs(X[order]) / N

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(t, np.real(x), lw=1)
    ax1.set(title=f"{title}\n(time domain)", xlabel="time (s)", ylabel="amplitude")
    ax1.grid(alpha=0.3)

    ax2.stem(f_plot, mag, basefmt=" ")
    ax2.set(title="magnitude spectrum\n(frequency domain)",
            xlabel="frequency (Hz)", ylabel="|X(f)|")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    # fig.savefig(fname, dpi=110)
    # plt.close(fig)
    plt.show()
    print(f"  saved {fname}")


# ----------------------------------------------------------------------------
# EXAMPLE 1: One pure sine wave -> one spike.
# The simplest possible fact: a single frequency lives in a single bin.
# This is the "Hello World" that makes the change-of-basis idea concrete.
# ----------------------------------------------------------------------------
def example_1_single_sine():
    print("Example 1: a single 5 Hz sine -> a single spike at 5 Hz")
    fs = 200          # sampling rate (Hz)
    T = 1.0           # total duration (s)  -> N = fs*T = 200 samples
    t = np.arange(0, T, 1/fs)
    f0 = 5            # signal frequency, chosen as an integer # of Hz on purpose
    x = np.sin(2*np.pi*f0*t)
    show_time_and_spectrum(t, x, fs, "5 Hz sine", "ex1_single_sine.png")
    # Takeaway: the spectrum has essentially one nonzero entry, at f0.


# ----------------------------------------------------------------------------
# EXAMPLE 2: A sum of sines -> a sum of spikes.
# LINEARITY: FFT(a+b) = FFT(a) + FFT(b). The transform reads a mixture as
# the literal list of its ingredients. This is *why* it is useful.
# ----------------------------------------------------------------------------
def example_2_sum_of_sines():
    print("Example 2: 3 Hz + 7 Hz + 12 Hz mixed together -> three spikes")
    fs = 200
    t = np.arange(0, 1.0, 1/fs)
    x = (1.0*np.sin(2*np.pi*3*t)      # amplitude 1.0 at 3 Hz
         + 0.5*np.sin(2*np.pi*7*t)    # amplitude 0.5 at 7 Hz
         + 0.8*np.sin(2*np.pi*12*t))  # amplitude 0.8 at 12 Hz
    show_time_and_spectrum(t, x, fs, "3+7+12 Hz mix", "ex2_sum_of_sines.png")
    # Takeaway: spike HEIGHTS recover the amplitudes (1.0, 0.5, 0.8).
    # The tangled time-domain wave is trivially separable in frequency.


# ----------------------------------------------------------------------------
# EXAMPLE 3: Phase. Magnitude alone is not the whole story.
# A sine and a cosine of the same frequency give the SAME magnitude spectrum
# but differ in PHASE (the angle of the complex number X[k]).
# Edge case for the "spectrum = frequency content" slogan: content includes
# phase, which the magnitude plot throws away.
# ----------------------------------------------------------------------------
def example_3_phase():
    print("Example 3: sine vs cosine -> same magnitude, different phase")
    fs = 200
    t = np.arange(0, 1.0, 1/fs)
    N = len(t)
    f0 = 5
    xs = np.sin(2*np.pi*f0*t)
    xc = np.cos(2*np.pi*f0*t)
    Xs, Xc = np.fft.fft(xs), np.fft.fft(xc)
    k = f0                                   # bin index for 5 Hz (since T=1s)
    print(f"  |Xs[{k}]| = {abs(Xs[k]):.1f},  |Xc[{k}]| = {abs(Xc[k]):.1f}  (equal magnitudes)")
    print(f"  angle(Xs[{k}]) = {np.degrees(np.angle(Xs[k])):7.1f} deg")
    print(f"  angle(Xc[{k}]) = {np.degrees(np.angle(Xc[k])):7.1f} deg  (differ by 90 deg)")
    # Takeaway: cos and sin are the same oscillation shifted by a quarter
    # period; that shift IS the phase. Real signal analysis needs both
    # magnitude and phase to be invertible.


# ----------------------------------------------------------------------------
# EXAMPLE 4: Spectral leakage.
# When the signal frequency is NOT an integer number of cycles in the window,
# energy spreads across many bins. This is the edge case that breaks the naive
# "one frequency -> one clean spike" reading. A window function tames it.
# ----------------------------------------------------------------------------
def example_4_leakage():
    print("Example 4: leakage when frequency falls between bins")
    fs = 200
    t = np.arange(0, 1.0, 1/fs)
    f0 = 5.5          # NOT an integer # of cycles in a 1 s window -> leaks
    x = np.sin(2*np.pi*f0*t)

    x_win = x * np.hanning(len(x))   # taper the ends to reduce leakage

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, sig, lab in [(axes[0], x, "no window"),
                         (axes[1], x_win, "Hann window")]:
        N = len(sig)
        f = np.fft.fftfreq(N, 1/fs)
        keep = f >= 0
        ax.stem(f[keep], np.abs(np.fft.fft(sig))[keep], basefmt=" ")
        ax.set(title=f"5.5 Hz, {lab}", xlabel="frequency (Hz)",
               ylabel="|X(f)|", xlim=(0, 20))
        ax.grid(alpha=0.3)
    fig.tight_layout()
    # fig.savefig("ex4_leakage.png", dpi=110)
    # plt.close(fig)
    plt.show()
    print("  saved ex4_leakage.png")
    # Takeaway: the DFT silently ASSUMES your window repeats periodically.
    # A 5.5-cycle window has a discontinuity at the seam -> broadband leakage.
    # Windowing trades a wider main lobe for much lower side lobes.


# ----------------------------------------------------------------------------
# EXAMPLE 5: Aliasing. Sampling too slowly makes a high frequency MASQUERADE
# as a low one. Nyquist theorem (a THEOREM, not a convention): to represent a
# signal whose highest frequency is B Hz, you must sample faster than 2B Hz.
# ----------------------------------------------------------------------------
def example_5_aliasing():
    print("Example 5: aliasing -- a 90 Hz tone sampled at 100 Hz looks like 10 Hz")
    fs = 100                       # Nyquist limit = 50 Hz
    f0 = 90                        # ABOVE Nyquist -> will alias
    t = np.arange(0, 1.0, 1/fs)
    x = np.sin(2*np.pi*f0*t)
    # densely-sampled "truth" for visual comparison
    t_fine = np.arange(0, 1.0, 1/2000)
    x_fine = np.sin(2*np.pi*f0*t_fine)
    alias = abs(((f0 + fs/2) % fs) - fs/2)   # predicted apparent frequency

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(t_fine, x_fine, alpha=0.3, label="true 90 Hz")
    ax1.plot(t, x, 'ro-', ms=3, lw=1, label=f"samples @ {fs} Hz")
    ax1.set(title="undersampling", xlabel="time (s)", xlim=(0, 0.2))
    ax1.legend(); ax1.grid(alpha=0.3)
    f = np.fft.fftfreq(len(x), 1/fs); keep = f >= 0
    ax2.stem(f[keep], np.abs(np.fft.fft(x))[keep], basefmt=" ")
    ax2.set(title=f"spectrum: peak appears at {alias:.0f} Hz, not 90",
            xlabel="frequency (Hz)", ylabel="|X(f)|")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    # fig.savefig("ex5_aliasing.png", dpi=110)
    # plt.close(fig)
    plt.show()
    print(f"  predicted alias frequency = {alias:.0f} Hz; saved ex5_aliasing.png")
    # Takeaway: frequencies above fs/2 fold back below it. This is why audio
    # gear low-pass filters BEFORE the analog-to-digital converter.


# ----------------------------------------------------------------------------
# EXAMPLE 6: A practical use -- denoising by zeroing small frequency bins.
# Pipeline: signal -> FFT -> kill bins below a threshold -> IFFT -> clean signal.
# Inverse transform IFFT undoes FFT exactly (up to floating point), so editing
# in the frequency domain and coming back is lossless except for what you cut.
# ----------------------------------------------------------------------------
def example_6_denoise():
    print("Example 6: denoise by thresholding in the frequency domain")
    fs = 500
    t = np.arange(0, 1.0, 1/fs)
    rng = np.random.default_rng(0)
    clean = np.sin(2*np.pi*8*t) + 0.5*np.sin(2*np.pi*20*t)
    noisy = clean + 1.2*rng.standard_normal(len(t))

    X = np.fft.fft(noisy)
    thresh = 0.1 * np.max(np.abs(X))   # keep only dominant bins
    X_clean = np.where(np.abs(X) > thresh, X, 0)
    recovered = np.real(np.fft.ifft(X_clean))

    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    for ax, sig, lab in [(axes[0], noisy, "noisy input"),
                         (axes[1], recovered, "recovered (FFT-thresholded)"),
                         (axes[2], clean, "true clean signal")]:
        ax.plot(t, sig, lw=1); ax.set(ylabel=lab); ax.grid(alpha=0.3)
    axes[-1].set_xlabel("time (s)")
    fig.suptitle("Denoising via the frequency domain")
    fig.tight_layout()
    # fig.savefig("ex6_denoise.png", dpi=110)
    # plt.close(fig)
    plt.show()
    err = np.sqrt(np.mean((recovered - clean)**2))
    print(f"  RMS error after denoising = {err:.3f}; saved ex6_denoise.png")
    # Takeaway: noise is spread thinly over all bins; real structure concentrates
    # in a few. Cutting low-magnitude bins keeps signal, drops noise.


if __name__ == "__main__":
    print("=" * 60)
    print("FOURIER TRANSFORM, ONE IDEA AT A TIME")
    print("=" * 60)
    example_1_single_sine()
    example_2_sum_of_sines()
    example_3_phase()
    example_4_leakage()
    example_5_aliasing()
    example_6_denoise()
    print("\nDone. Open the PNGs in order ex1..ex6 to follow the logic chain.")
