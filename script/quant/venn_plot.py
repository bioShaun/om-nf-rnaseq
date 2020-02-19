#! /usr/bin/env python

import os
import sys
import math
import click
import pandas as pd
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
from itertools import chain
from collections import Iterable
mpl.use('Agg')


DIFF_LIST_SFX = 'edgeR.DE_results.diffgenes.txt'


default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [255, 140, 0, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


def save_mkdir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def write_obj_to_file(obj, fn, append=False):
    fh = open(fn, 'a' if append is True else 'w')
    if type(obj) is str:
        fh.write('%s\n' % obj)
    elif type(obj) is list or type(obj) is set:
        for item in obj:
            fh.write('%s\n' % item)
    elif type(obj) is dict:
        for key, val in obj.iteritems():
            fh.write('%s\t%s\n' % (key, val))
    else:
        raise TypeError('invalid type for %s' % obj)
    fh.close()


def draw_ellipse(fig, ax, x, y, w, h, a, fillcolor):
    e = patches.Ellipse(
        xy=(x, y),
        width=w,
        height=h,
        angle=a,
        color=fillcolor)
    ax.add_patch(e)


def draw_triangle(fig, ax, x1, y1, x2, y2, x3, y3, fillcolor):
    xy = [
        (x1, y1),
        (x2, y2),
        (x3, y3),
    ]
    polygon = patches.Polygon(
        xy=xy,
        closed=True,
        color=fillcolor)
    ax.add_patch(polygon)


def draw_text(fig, ax, x, y, text, color=[0, 0, 0, 1]):
    ax.text(
        x, y, text,
        horizontalalignment='center',
        verticalalignment='center',
        fontsize=14,
        color=color)


def draw_annotate(fig, ax, x, y, textx, texty, text, color=[0, 0, 0, 1], arrowcolor=[0, 0, 0, 0.3]):
    plt.annotate(
        text,
        xy=(x, y),
        xytext=(textx, texty),
        arrowprops=dict(color=arrowcolor, shrink=0, width=0.5, headwidth=8),
        fontsize=14,
        color=color,
        xycoords="data",
        textcoords="data",
        horizontalalignment='center',
        verticalalignment='center'
    )


def get_labels(data, fill=["number"]):
    """
    get a dict of labels for groups in data

    @type data: list[Iterable]
    @rtype: dict[str, str]

    input
      data: data to get label for
      fill: ["number"|"logic"|"percent"]

    return
      labels: a dict of labels for different sets

    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill=["number"])
    Out[12]:
    {'001': '0',
     '010': '5',
     '011': '0',
     '100': '3',
     '101': '2',
     '110': '2',
     '111': '3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i]
                                 for i in range(N) if key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if "logic" in fill:
        for k in set_collections:
            labels[k] = k + ": "
    if "number" in fill:
        for k in set_collections:
            labels[k] += str(len(set_collections[k]))
    if "percent" in fill:
        data_size = len(s_all)
        for k in set_collections:
            labels[k] += "(%.1f%%)" % (100.0 *
                                       len(set_collections[k]) / data_size)

    return labels, set_collections


def venn2(labels, names=['A', 'B'], **options):
    """
    plots a 2-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('01', '10', '11'),
              hence a valid set could look like: {'01': 'text 1', '10': 'text 2', '11': 'text 3'}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(2)])
    figsize = options.get('figsize', (9, 7))
    dpi = options.get('dpi', 96)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=0.7)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.375, 0.3, 0.5, 0.5, 0.0, colors[0])
    draw_ellipse(fig, ax, 0.625, 0.3, 0.5, 0.5, 0.0, colors[1])
    draw_text(fig, ax, 0.74, 0.30, labels.get('01', ''))
    draw_text(fig, ax, 0.26, 0.30, labels.get('10', ''))
    draw_text(fig, ax, 0.50, 0.30, labels.get('11', ''))

    # legend
    draw_text(fig, ax, 0.20, 0.56, names[0], colors[0])
    draw_text(fig, ax, 0.80, 0.56, names[1], colors[1])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn3(labels, names=['A', 'B', 'C'], **options):
    """
    plots a 3-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('001', '010', '100', ...),
              hence a valid set could look like: {'001': 'text 1', '010': 'text 2', '100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(3)])
    figsize = options.get('figsize', (9, 9))
    dpi = options.get('dpi', 96)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.333, 0.633, 0.5, 0.5, 0.0, colors[0])
    draw_ellipse(fig, ax, 0.666, 0.633, 0.5, 0.5, 0.0, colors[1])
    draw_ellipse(fig, ax, 0.500, 0.310, 0.5, 0.5, 0.0, colors[2])
    draw_text(fig, ax, 0.50, 0.27, labels.get('001', ''))
    draw_text(fig, ax, 0.73, 0.65, labels.get('010', ''))
    draw_text(fig, ax, 0.61, 0.46, labels.get('011', ''))
    draw_text(fig, ax, 0.27, 0.65, labels.get('100', ''))
    draw_text(fig, ax, 0.39, 0.46, labels.get('101', ''))
    draw_text(fig, ax, 0.50, 0.65, labels.get('110', ''))
    draw_text(fig, ax, 0.50, 0.51, labels.get('111', ''))

    # legend
    draw_text(fig, ax, 0.15, 0.87, names[0], colors[0])
    draw_text(fig, ax, 0.85, 0.87, names[1], colors[1])
    draw_text(fig, ax, 0.50, 0.02, names[2], colors[2])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn4(labels, names=['A', 'B', 'C', 'D'], **options):
    """
    plots a 4-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('0001', '0010', '0100', ...),
              hence a valid set could look like: {'0001': 'text 1', '0010': 'text 2', '0100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(4)])
    figsize = options.get('figsize', (12, 12))
    dpi = options.get('dpi', 96)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.350, 0.400, 0.72, 0.45, 140.0, colors[0])
    draw_ellipse(fig, ax, 0.450, 0.500, 0.72, 0.45, 140.0, colors[1])
    draw_ellipse(fig, ax, 0.544, 0.500, 0.72, 0.45, 40.0, colors[2])
    draw_ellipse(fig, ax, 0.644, 0.400, 0.72, 0.45, 40.0, colors[3])
    draw_text(fig, ax, 0.85, 0.42, labels.get('0001', ''))
    draw_text(fig, ax, 0.68, 0.72, labels.get('0010', ''))
    draw_text(fig, ax, 0.77, 0.59, labels.get('0011', ''))
    draw_text(fig, ax, 0.32, 0.72, labels.get('0100', ''))
    draw_text(fig, ax, 0.71, 0.30, labels.get('0101', ''))
    draw_text(fig, ax, 0.50, 0.66, labels.get('0110', ''))
    draw_text(fig, ax, 0.65, 0.50, labels.get('0111', ''))
    draw_text(fig, ax, 0.14, 0.42, labels.get('1000', ''))
    draw_text(fig, ax, 0.50, 0.17, labels.get('1001', ''))
    draw_text(fig, ax, 0.29, 0.30, labels.get('1010', ''))
    draw_text(fig, ax, 0.39, 0.24, labels.get('1011', ''))
    draw_text(fig, ax, 0.23, 0.59, labels.get('1100', ''))
    draw_text(fig, ax, 0.61, 0.24, labels.get('1101', ''))
    draw_text(fig, ax, 0.35, 0.50, labels.get('1110', ''))
    draw_text(fig, ax, 0.50, 0.38, labels.get('1111', ''))

    # legend
    draw_text(fig, ax, 0.13, 0.18, names[0], colors[0])
    draw_text(fig, ax, 0.18, 0.83, names[1], colors[1])
    draw_text(fig, ax, 0.82, 0.83, names[2], colors[2])
    draw_text(fig, ax, 0.87, 0.18, names[3], colors[3])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn5(labels, names=['A', 'B', 'C', 'D', 'E'], **options):
    """
    plots a 5-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('00001', '00010', '00100', ...),
              hence a valid set could look like: {'00001': 'text 1', '00010': 'text 2', '00100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(5)])
    figsize = options.get('figsize', (13, 13))
    dpi = options.get('dpi', 96)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.0, top=1.0)
    ax.set_xlim(left=0.0, right=1.0)

    # body
    draw_ellipse(fig, ax, 0.428, 0.449, 0.87, 0.50, 155.0, colors[0])
    draw_ellipse(fig, ax, 0.469, 0.543, 0.87, 0.50, 82.0, colors[1])
    draw_ellipse(fig, ax, 0.558, 0.523, 0.87, 0.50, 10.0, colors[2])
    draw_ellipse(fig, ax, 0.578, 0.432, 0.87, 0.50, 118.0, colors[3])
    draw_ellipse(fig, ax, 0.489, 0.383, 0.87, 0.50, 46.0, colors[4])
    draw_text(fig, ax, 0.27, 0.11, labels.get('00001', ''))
    draw_text(fig, ax, 0.72, 0.11, labels.get('00010', ''))
    draw_text(fig, ax, 0.55, 0.13, labels.get('00011', ''))
    draw_text(fig, ax, 0.91, 0.58, labels.get('00100', ''))
    draw_text(fig, ax, 0.78, 0.64, labels.get('00101', ''))
    draw_text(fig, ax, 0.84, 0.41, labels.get('00110', ''))
    draw_text(fig, ax, 0.76, 0.55, labels.get('00111', ''))
    draw_text(fig, ax, 0.51, 0.90, labels.get('01000', ''))
    draw_text(fig, ax, 0.39, 0.15, labels.get('01001', ''))
    draw_text(fig, ax, 0.42, 0.78, labels.get('01010', ''))
    draw_text(fig, ax, 0.50, 0.15, labels.get('01011', ''))
    draw_text(fig, ax, 0.67, 0.76, labels.get('01100', ''))
    draw_text(fig, ax, 0.70, 0.71, labels.get('01101', ''))
    draw_text(fig, ax, 0.51, 0.74, labels.get('01110', ''))
    draw_text(fig, ax, 0.64, 0.67, labels.get('01111', ''))
    draw_text(fig, ax, 0.10, 0.61, labels.get('10000', ''))
    draw_text(fig, ax, 0.20, 0.31, labels.get('10001', ''))
    draw_text(fig, ax, 0.76, 0.25, labels.get('10010', ''))
    draw_text(fig, ax, 0.65, 0.23, labels.get('10011', ''))
    draw_text(fig, ax, 0.18, 0.50, labels.get('10100', ''))
    draw_text(fig, ax, 0.21, 0.37, labels.get('10101', ''))
    draw_text(fig, ax, 0.81, 0.37, labels.get('10110', ''))
    draw_text(fig, ax, 0.74, 0.40, labels.get('10111', ''))
    draw_text(fig, ax, 0.27, 0.70, labels.get('11000', ''))
    draw_text(fig, ax, 0.34, 0.25, labels.get('11001', ''))
    draw_text(fig, ax, 0.33, 0.72, labels.get('11010', ''))
    draw_text(fig, ax, 0.51, 0.22, labels.get('11011', ''))
    draw_text(fig, ax, 0.25, 0.58, labels.get('11100', ''))
    draw_text(fig, ax, 0.28, 0.39, labels.get('11101', ''))
    draw_text(fig, ax, 0.36, 0.66, labels.get('11110', ''))
    draw_text(fig, ax, 0.51, 0.47, labels.get('11111', ''))

    # legend
    draw_text(fig, ax, 0.02, 0.72, names[0], colors[0])
    draw_text(fig, ax, 0.72, 0.94, names[1], colors[1])
    draw_text(fig, ax, 0.97, 0.74, names[2], colors[2])
    draw_text(fig, ax, 0.88, 0.05, names[3], colors[3])
    draw_text(fig, ax, 0.12, 0.05, names[4], colors[4])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def venn6(labels, names=['A', 'B', 'C', 'D', 'E'], **options):
    """
    plots a 6-set Venn diagram

    @type labels: dict[str, str]
    @type names: list[str]
    @rtype: (Figure, AxesSubplot)

    input
      labels: a label dict where keys are identified via binary codes ('000001', '000010', '000100', ...),
              hence a valid set could look like: {'000001': 'text 1', '000010': 'text 2', '000100': 'text 3', ...}.
              unmentioned codes are considered as ''.
      names:  group names
      more:   colors, figsize, dpi

    return
      pyplot Figure and AxesSubplot object
    """
    colors = options.get('colors', [default_colors[i] for i in range(6)])
    figsize = options.get('figsize', (20, 20))
    dpi = options.get('dpi', 96)

    fig = plt.figure(0, figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    ax.set_axis_off()
    ax.set_ylim(bottom=0.230, top=0.845)
    ax.set_xlim(left=0.173, right=0.788)

    # body
    draw_triangle(fig, ax, 0.637, 0.921, 0.649, 0.274, 0.188, 0.667, colors[0])
    draw_triangle(fig, ax, 0.981, 0.769, 0.335, 0.191, 0.393, 0.671, colors[1])
    draw_triangle(fig, ax, 0.941, 0.397, 0.292, 0.475, 0.456, 0.747, colors[2])
    draw_triangle(fig, ax, 0.662, 0.119, 0.316, 0.548, 0.662, 0.700, colors[3])
    draw_triangle(fig, ax, 0.309, 0.081, 0.374, 0.718, 0.681, 0.488, colors[4])
    draw_triangle(fig, ax, 0.016, 0.626, 0.726, 0.687, 0.522, 0.327, colors[5])
    draw_text(fig, ax, 0.212, 0.562, labels.get('000001', ''))
    draw_text(fig, ax, 0.430, 0.249, labels.get('000010', ''))
    draw_text(fig, ax, 0.356, 0.444, labels.get('000011', ''))
    draw_text(fig, ax, 0.609, 0.255, labels.get('000100', ''))
    draw_text(fig, ax, 0.323, 0.546, labels.get('000101', ''))
    draw_text(fig, ax, 0.513, 0.316, labels.get('000110', ''))
    draw_text(fig, ax, 0.523, 0.348, labels.get('000111', ''))
    draw_text(fig, ax, 0.747, 0.458, labels.get('001000', ''))
    draw_text(fig, ax, 0.325, 0.492, labels.get('001001', ''))
    draw_text(fig, ax, 0.670, 0.481, labels.get('001010', ''))
    draw_text(fig, ax, 0.359, 0.478, labels.get('001011', ''))
    draw_text(fig, ax, 0.653, 0.444, labels.get('001100', ''))
    draw_text(fig, ax, 0.344, 0.526, labels.get('001101', ''))
    draw_text(fig, ax, 0.653, 0.466, labels.get('001110', ''))
    draw_text(fig, ax, 0.363, 0.503, labels.get('001111', ''))
    draw_text(fig, ax, 0.750, 0.616, labels.get('010000', ''))
    draw_text(fig, ax, 0.682, 0.654, labels.get('010001', ''))
    draw_text(fig, ax, 0.402, 0.310, labels.get('010010', ''))
    draw_text(fig, ax, 0.392, 0.421, labels.get('010011', ''))
    draw_text(fig, ax, 0.653, 0.691, labels.get('010100', ''))
    draw_text(fig, ax, 0.651, 0.644, labels.get('010101', ''))
    draw_text(fig, ax, 0.490, 0.340, labels.get('010110', ''))
    draw_text(fig, ax, 0.468, 0.399, labels.get('010111', ''))
    draw_text(fig, ax, 0.692, 0.545, labels.get('011000', ''))
    draw_text(fig, ax, 0.666, 0.592, labels.get('011001', ''))
    draw_text(fig, ax, 0.665, 0.496, labels.get('011010', ''))
    draw_text(fig, ax, 0.374, 0.470, labels.get('011011', ''))
    draw_text(fig, ax, 0.653, 0.537, labels.get('011100', ''))
    draw_text(fig, ax, 0.652, 0.579, labels.get('011101', ''))
    draw_text(fig, ax, 0.653, 0.488, labels.get('011110', ''))
    draw_text(fig, ax, 0.389, 0.486, labels.get('011111', ''))
    draw_text(fig, ax, 0.553, 0.806, labels.get('100000', ''))
    draw_text(fig, ax, 0.313, 0.604, labels.get('100001', ''))
    draw_text(fig, ax, 0.388, 0.694, labels.get('100010', ''))
    draw_text(fig, ax, 0.375, 0.633, labels.get('100011', ''))
    draw_text(fig, ax, 0.605, 0.359, labels.get('100100', ''))
    draw_text(fig, ax, 0.334, 0.555, labels.get('100101', ''))
    draw_text(fig, ax, 0.582, 0.397, labels.get('100110', ''))
    draw_text(fig, ax, 0.542, 0.372, labels.get('100111', ''))
    draw_text(fig, ax, 0.468, 0.708, labels.get('101000', ''))
    draw_text(fig, ax, 0.355, 0.572, labels.get('101001', ''))
    draw_text(fig, ax, 0.420, 0.679, labels.get('101010', ''))
    draw_text(fig, ax, 0.375, 0.597, labels.get('101011', ''))
    draw_text(fig, ax, 0.641, 0.436, labels.get('101100', ''))
    draw_text(fig, ax, 0.348, 0.538, labels.get('101101', ''))
    draw_text(fig, ax, 0.635, 0.453, labels.get('101110', ''))
    draw_text(fig, ax, 0.370, 0.548, labels.get('101111', ''))
    draw_text(fig, ax, 0.594, 0.689, labels.get('110000', ''))
    draw_text(fig, ax, 0.579, 0.670, labels.get('110001', ''))
    draw_text(fig, ax, 0.398, 0.670, labels.get('110010', ''))
    draw_text(fig, ax, 0.395, 0.653, labels.get('110011', ''))
    draw_text(fig, ax, 0.633, 0.682, labels.get('110100', ''))
    draw_text(fig, ax, 0.616, 0.656, labels.get('110101', ''))
    draw_text(fig, ax, 0.587, 0.427, labels.get('110110', ''))
    draw_text(fig, ax, 0.526, 0.415, labels.get('110111', ''))
    draw_text(fig, ax, 0.495, 0.677, labels.get('111000', ''))
    draw_text(fig, ax, 0.505, 0.648, labels.get('111001', ''))
    draw_text(fig, ax, 0.428, 0.663, labels.get('111010', ''))
    draw_text(fig, ax, 0.430, 0.631, labels.get('111011', ''))
    draw_text(fig, ax, 0.639, 0.524, labels.get('111100', ''))
    draw_text(fig, ax, 0.591, 0.604, labels.get('111101', ''))
    draw_text(fig, ax, 0.622, 0.477, labels.get('111110', ''))
    draw_text(fig, ax, 0.501, 0.523, labels.get('111111', ''))

    # legend
    draw_text(fig, ax, 0.674, 0.824, names[0], colors[0])
    draw_text(fig, ax, 0.747, 0.751, names[1], colors[1])
    draw_text(fig, ax, 0.739, 0.396, names[2], colors[2])
    draw_text(fig, ax, 0.700, 0.247, names[3], colors[3])
    draw_text(fig, ax, 0.291, 0.255, names[4], colors[4])
    draw_text(fig, ax, 0.203, 0.484, names[5], colors[5])
    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    return fig, ax


def output_plot(fig, out_prefix, pdf=True):
    png_file = '{p}.png'.format(p=out_prefix)
    fig.savefig(png_file, dpi=300)
    if pdf:
        pdf_file = '{p}.pdf'.format(p=out_prefix)
        fig.savefig(pdf_file)
    fig.clear()
    plt.close(fig)


VENN_FUNCS = {key: val for key, val in locals().items() if 'venn' in key}


def plot_one_combination(diff_dir, combination, out_dir):
    combination_list = combination.split(',')
    comb_num = len(combination_list)
    # check combination number, only support 2-5 combination to draw venn
    if comb_num < 2 or comb_num > 5:
        click.echo('Sets number error.\nSupported 2-5 sets venn plot.')
        sys.exit(1)

    # read diff-gene list
    gene_list = []
    for each_com in combination_list:
        each_com_diff_genes = os.path.join(
            diff_dir, each_com, '{c}.ALL.{s}'.format(
                c=each_com, s=DIFF_LIST_SFX))
        diff_df = pd.read_table(each_com_diff_genes, header=None, index_col=0)
        gene_list.append(diff_df.index)

    # draw venn, output 2 figs, one with logical numbers and the other without
    labels, col_set = get_labels(gene_list, fill=['number', 'logic'])
    labels_simple, col_set = get_labels(gene_list, fill=['number'])
    venn_prefix = '__'.join(combination_list)
    result_plt_dir = os.path.join(
        out_dir, 'venn_plot', venn_prefix)
    save_mkdir(result_plt_dir)
    report_plt_dir = os.path.join(
        out_dir, 'report_plot',
    )
    save_mkdir(report_plt_dir)
    venn_plot = VENN_FUNCS[f'venn{comb_num}']
    fig_detail, ax = venn_plot(labels, names=combination_list)
    venn_detail = os.path.join(result_plt_dir, 'venn_plot.detail')
    output_plot(fig_detail, venn_detail)
    fig_simple, ax = venn_plot(labels_simple, names=combination_list)
    venn_simple = os.path.join(result_plt_dir, 'venn_plot')
    output_plot(fig_simple, venn_simple)
    fig_simple2, ax = venn_plot(labels_simple, names=combination_list)
    venn_report = os.path.join(report_plt_dir, f'{venn_prefix}.venn')
    output_plot(fig_simple2, venn_report, pdf=False)

    # output gene ids in each part of venn
    for each_part in col_set:
        out_file = os.path.join(result_plt_dir, '{n}.txt'.format(n=each_part))
        write_obj_to_file(col_set[each_part], out_file)


@click.command()
@click.option('-d', '--diff_dir', type=click.Path(exists=True), required=True,
              help='differential analysis directory.')
@click.option('-s', '--combination', type=click.STRING,
              help='venn plot sets, limit is 2 ~ 5 Sets, seperated with ",".')
@click.option('-f', '--combination_file', type=click.Path(exists=True),
              help='venn plot sets file')
@click.option('-o', '--out_dir', type=click.Path(), required=True,
              help='output directory.')
@click.option('-a', '--all_combine', is_flag=True,
              help='run all possible combination between 2-5 sets of sets \
              provided.')
def main(diff_dir, combination_file, combination, out_dir, all_combine):
    if combination_file:
        combinations = [each.strip() for each in open(combination_file)]
        [plot_one_combination(diff_dir, cmb, out_dir) for cmb in combinations]
    elif combination:
        if not all_combine:
            plot_one_combination(diff_dir, combination, out_dir)
        else:
            combination_list = combination.split(',')
            max_num = min(6, len(combination_list) + 1)
            for each_num in range(2, max_num):
                all_com = itertools.combinations(combination_list, each_num)
                for each_com in all_com:
                    each_com_name = ','.join(each_com)
                    plot_one_combination(diff_dir, each_com_name, out_dir)
    else:
        click.echo('--combination_file or --combination is needed')
        sys.exit(1)


if __name__ == '__main__':
    main()
