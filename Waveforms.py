__author__ = 'Martijn Schut'

from psychopy import sound
import math

wave = sound.Sound(value = 500)

def sine_wave(freq = 100, sample_rate = 44100, duration = 1):
    """

    Produces a trianglewave PsychoPy/PyGame object. Linearly interpolates between
    the minimum and the maximum.

    freq: frequency in Hz
    sample_rate: well, sample rate
    duration: duration of the sound in seconds

    returns: a psychopy sound object
    """

    duration_in_samples = duration * sample_rate
    Hz = [math.sin(2.0 * math.pi * freq * tp / sample_rate) for tp in xrange(0, duration_in_samples)]

    wave = sound.SoundPygame(value = Hz)

    return wave
