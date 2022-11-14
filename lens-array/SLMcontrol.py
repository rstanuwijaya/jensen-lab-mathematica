from math import ceil, floor, pi, acos, prod, sqrt, sin, cos, tan
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import cv2

# import detect_heds_module_path
# from holoeye import slmdisplaysdk

# slm = slmdisplaysdk.SLMInstance()

# # Check if the library implements the required version
# if not slm.requiresVersion(3):
#     exit(1)


from colorsys import hls_to_rgb
def colorize(z):
    n,m = z.shape
    c = np.zeros((n,m,3))
    c[np.isinf(z)] = (1.0, 1.0, 1.0)
    c[np.isnan(z)] = (0.5, 0.5, 0.5)

    idx = ~(np.isinf(z) + np.isnan(z))
    A = (np.angle(z[idx]) + np.pi) / (2*np.pi)
    A = (A + 0.5) % 1.0
    B = 1.0 - 1.0/(1.0+abs(z[idx])**0.3)
    c[idx] = [hls_to_rgb(a, b, 0.8) for a,b in zip(A,B)]
    return c

slm_res = 8
slm_dim = np.array((1080, 1920))
lens_size = 64
wl = 0.810
deflect_angle = 0.5
offset = np.array([0, 0])
pattern_size, redun_size = np.divmod(slm_dim, lens_size)
k0 = 2*pi/wl

def slm_phaseamp_control(profile) -> np.ndarray:
    phase_plus = np.angle(profile) + np.arccos(np.abs(profile))
    phase_minus = np.angle(profile) - np.arccos(np.abs(profile))
    meta = np.array([[np.exp(1j*phase_plus), np.exp(1j*phase_minus)],
                     [np.exp(1j*phase_minus), np.exp(1j*phase_plus)]])
    meta = np.transpose(meta, (2, 0, 3, 1))
    meta = np.repeat(meta, lens_size//2, axis=0)
    meta = np.repeat(meta, lens_size//2, axis=2)
    meta = np.reshape(
        meta, (meta.shape[0]*meta.shape[1], meta.shape[2]*meta.shape[3]))
    return meta

def slm_constant_control(profile) -> np.ndarray:
    meta = np.array(np.exp(1j*np.angle(profile)))
    meta = np.repeat(meta, lens_size, axis=0)
    meta = np.repeat(meta, lens_size, axis=1)
    # meta = np.reshape(
    #     meta, (meta.shape[0]*meta.shape[1], meta.shape[2]*meta.shape[3]))
    return meta

def create_phase_gradient() -> np.ndarray:
    x = np.arange(0, slm_dim[0]*slm_res, slm_res)
    y = np.arange(0, slm_dim[1]*slm_res, slm_res)
    y, x = np.meshgrid(y, x)
    pg = np.exp(1j*k0*sin(deflect_angle)*x)
    return pg


def slm_gradient_control(profile) -> np.ndarray:
    constant = np.ones(
        (pattern_size[0]*lens_size, pattern_size[1]*lens_size), dtype=complex)
    gradient = phase_gradient[:pattern_size[0]*lens_size, :pattern_size[1]*lens_size]
    cond = np.kron(profile, np.ones((lens_size, lens_size)))
    meta = np.where(cond, gradient, constant)
    return meta


phase_gradient = create_phase_gradient()


def create_pattern(pattern, offset0, offset1, lens_size0, deflect_angle0):
    global profile, offset, lens_size, pattern_size, redun_size, deflect_angle, phase_gradient
    profile = np.zeros(pattern_size, dtype=complex)
    pattern = np.array(pattern)
    profile[offset[0]:offset[0]+pattern.shape[0],
            offset[1]:offset[1]+pattern.shape[1]] = pattern

    offset[0] = int(offset0)
    offset[1] = int(offset1)
    lens_size = lens_size0
    pattern_size, redun_size = np.divmod(slm_dim, lens_size)
    if deflect_angle != deflect_angle0:
        deflect_angle = deflect_angle0
        phase_gradient = create_phase_gradient()

    profile = np.zeros(pattern_size, dtype=complex)
    profile[offset[0]: offset[0]+pattern.shape[0],
            offset[1]: offset[1]+pattern.shape[1]] = pattern
    try:
        slm_phase = np.angle(slm_constant_control(
            profile)*slm_gradient_control(np.abs(profile) == 0))
        slm_phase = np.pad(
            slm_phase, ((0, redun_size[0]), (0, redun_size[1])), 'constant')
    except Exception as e:
        slm_phase = np.zeros(slm_dim)
        # raise (e)
    return slm_phase


def display_slm(slm_phase):
    return True