from psychopy import visual
from psychopy import sound
from psychopy import core
import math
from scipy import signal
import numpy as np
try:
    from pyfftw.interfaces.numpy_fft import rfft, irfft       # Performs much better than numpy's fftpack
except ImportError:                                    # Use monkey-patching np.fft perhaps instead?
    from numpy.fft import rfft, irfft

def ms(x):
    """Mean value of signal `x` squared.
    :param x: Dynamic quantity.
    :returns: Mean squared of `x`.
    """
    return (np.abs(x)**2.0).mean()

def normalize(y, x=None):
    """normalize power in y to a (standard normal) white noise signal.
    Optionally normalize to power in signal `x`.
    #The mean power of a Gaussian with :math:`\\mu=0` and :math:`\\sigma=1` is 1.
    """
    #return y * np.sqrt( (np.abs(x)**2.0).mean() / (np.abs(y)**2.0).mean() )
    if x is not None:
        x = ms(x)
    else:
        x = 1.0
    return y * np.sqrt( x / ms(y) )

def noise_spectrum(N, spectrum = 1):
    """ Spectrum indicates by what you times the spectrum, -1 generates pink noise
    """
    uneven = N % 2
    X = np.random.randn(N // 2 + 1 + uneven) + 1j * np.random.randn(N // 2 + 1 + uneven)
    S = np.sqrt(np.arange(len(X)) + 1.)  # +1 to avoid divide by zero
    y = (irfft((X * spectrum) * S)).real
    if uneven:
        y = y[:-1]
    return normalize(y)

def pink(N):
    """
    Pink noise.

    :param N: Amount of samples.

    Pink noise has equal power in bands that are proportionally wide.
    Power density decreases with 3 dB per octave.

    """
    # This method uses the filter with the following coefficients.
    # b = np.array([0.049922035, -0.095993537, 0.050612699, -0.004408786])
    # a = np.array([1, -2.494956002, 2.017265875, -0.522189400])
    # return lfilter(B, A, np.random.randn(N))
    # Another way would be using the FFT
    # x = np.random.randn(N)
    # X = rfft(x) / N
    uneven = N % 2
    X = np.random.randn(N // 2 + 1 + uneven) + 1j * np.random.randn(N // 2 + 1 + uneven)
    S = np.sqrt(np.arange(len(X)) + 1.)  # +1 to avoid divide by zero
    y = (irfft(X / S)).real
    if uneven:
        y = y[:-1]
    return normalize(y)


def brown(N):
    """
    Violet noise.

    :param N: Amount of samples.

    Power decreases with -3 dB per octave.
    Power density decreases with 6 dB per octave.
    """
    uneven = N % 2
    X = np.random.randn(N // 2 + 1 + uneven) + 1j * np.random.randn(N // 2 + 1 + uneven)
    S = (np.arange(len(X)))  # Filter
    y = (irfft(X * S)).real
    if uneven:
        y = y[:-1]
    return normalize(y)

def blue(N):
    """
    Blue noise.

    :param N: Amount of samples.

    Power increases with 6 dB per octave.
    Power density increases with 3 dB per octave.

    """
    uneven = N % 2
    X = np.random.randn(N // 2 + 1 + uneven) + 1j * np.random.randn(N // 2 + 1 + uneven)
    S = np.sqrt(np.arange(len(X)))  # Filter
    y = (irfft(X * S)).real
    if uneven:
        y = y[:-1]
    return normalize(y)

class waveform(object):
    def __init__(self, wavetype = 'sine', wave_spectrum = None, freq = 1600, sample_rate = 44100, duration = 1.0, pan = float(0.0), ITD = False, weighted = True):
        self.freq           = freq
        self.sample_rate    = sample_rate
        self.duration       = float(duration)
        self.pan            = float(pan)
        self.weight         = 1
        self.weighted       = weighted
        self.wavetype       = wavetype
        self.rampdur        = 3000
        self.ramp_up        = np.logspace(-3, 0, num = self.rampdur, endpoint=True)
        self.ramp_down      = np.logspace(0, -3, num = self.rampdur, endpoint=True)
        for x in range(0, len(self.ramp_down)):
            if self.ramp_down[x] < 0:
                self.ramp_down[x] = 0
        for x in range(0, len(self.ramp_up)):
            if self.ramp_up[x] < 0:
                self.ramp_up[x] = 0
        duration_in_samples = duration * sample_rate

        if self.pan > 1.0:
            self.pan = 1.0
        if self.pan < -1.0:
            self.pan = -1.0

        Lw                  = (0.5 + (-self.pan / 2))
        Rw                  = (0.5 + (self.pan / 2))

        if Lw > 1.0:
            Lw = 1.0
        elif Lw < 0.0:
            Lw = 0
        if Rw > 1.0:
            Rw = 1.0
        elif Rw < 0.0:
            Rw = 0.0
        print Rw, Lw

        from numpy import zeros
        self.stereoHz = zeros((duration_in_samples, 2))

        t = np.linspace(0, 1, duration * sample_rate, endpoint=False)
        if wavetype == 'sine' and wave_spectrum is None:

            self.stereoHz[0:self.rampdur - 1, 0] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Lw * self.weight * self.ramp_up[tp] for tp in xrange(0, 2999)]
            self.stereoHz[0:self.rampdur - 1, 1] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Rw * self.weight * self.ramp_up[tp] for tp in xrange(0, 2999)]

            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur , 0] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Lw * self.weight for tp in xrange(3000, int(round(duration_in_samples)) - 3000)]
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur , 1] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Rw * self.weight for tp in xrange(3000, int(round(duration_in_samples)) - 3000)]

            self.stereoHz[duration_in_samples - self.rampdur:duration_in_samples, 0] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Lw * self.weight * self.ramp_down[tp] for tp in xrange(0, 3000)]
            self.stereoHz[duration_in_samples - self.rampdur:duration_in_samples, 1] = [math.sin(2.0 * math.pi * self.freq * tp / self.sample_rate) * Rw * self.weight * self.ramp_down[tp] for tp in xrange(0, 3000)]

        if wavetype == 'wn' and wave_spectrum is None:
            s = np.random.sample(size = duration_in_samples)
            s *= 2
            s -= 1
            self.stereoHz[0:duration_in_samples, 0] = [sample * Lw for sample in s]
            self.stereoHz[0:duration_in_samples, 1] = [sample * Rw for sample in s]

        if wavetype == 'pn' and wave_spectrum is None:
            s = pink(duration_in_samples)
            for i in range(0, len(s)):
                if s[i] < -1:
                    s[i] = -1
                if s[i] > 1:
                    s[i] = 1

            self.stereoHz[0:duration_in_samples, 0] = [sample * Lw for sample in s]
            self.stereoHz[0:duration_in_samples, 1] = [sample * Rw for sample in s]

        if wavetype == 'bn' and wave_spectrum is None:
            s = brown(duration_in_samples)
            for i in range(0, len(s)):
                if s[i] < -1:
                    s[i] = -1
                if s[i] > 1:
                    s[i] = 1
            self.stereoHz[0:duration_in_samples, 0] = [sample * Lw for sample in s]
            self.stereoHz[0:duration_in_samples, 1] = [sample * Rw for sample in s]

        elif wavetype == 'square' and wave_spectrum is None:
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 0] = signal.square(2.0 * np.pi * self.freq * t) * Lw * self.weight * 0.5
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 1] = signal.square(2.0 * np.pi * self.freq * t) * Rw * self.weight * 0.5
        elif wavetype == 'sawtooth' and wave_spectrum is None:
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 0] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 1) * Lw * self.weight * 0.5
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 1] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 1) * Rw * self.weight * 0.5
        elif wavetype == 'sawtoothd' and wave_spectrum is None:
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 0] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 0) * Lw * self.weight * 0.5
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 1] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 0) * Rw * self.weight * 0.5
        elif wavetype == 'triangle' and wave_spectrum is None:
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 0] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 0.5) * Lw * self.weight * 0.5
            self.stereoHz[self.rampdur:duration_in_samples - self.rampdur, 1] = signal.sawtooth(2.0 * np.pi * self.freq * t, width = 0.5) * Rw * self.weight * 0.5
        elif wave_spectrum is not None:
            s = noise_spectrum(duration_in_samples, spectrum=wave_spectrum)
            for i in range(0, len(s)):
                if s[i] < -1:
                    s[i] = -1
                if s[i] > 1:
                    s[i] = 1

            self.stereoHz[0:duration_in_samples, 0] = [sample * Lw for sample in s]
            self.stereoHz[0:duration_in_samples, 1] = [sample * Rw for sample in s]
        else:
            raise Exception
        self.wave = sound.SoundPygame(value = self.stereoHz)

    def play(self):
        self.wave.play()

    def stop(self):
        self.wave.stop()

    def get_pygame_sound(self, Hz = None):
        if Hz is None:
            wave = sound.SoundPygame(value = self.Hz) 
        else:
            wave =  sound.SoundPygame(value = Hz)

        return wave

