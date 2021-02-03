---
layout: post
title:  "BigDFT Tutorials"
date:   2021-02-01 01:00:00 +0900
categories: blog
draft: false
---

Recently, we have been creating a new version of the website for the BigDFT
code. One of the goals for the new website is to present a new workflow
for using BigDFT. In the past, the code was used by preparing input files
manually and then running the program. However, for a few years now we
have been developing the PyBigDFT python bindings for BigDFT [1]. The new
workflow will have you write python code that can handle everything from 
pre-processing, to running calculations, to analyzing the results.

With this in mind, we have been developing some new tutorials for how to
use the code. Each of these tutorials is build from a jupyter notebook that
we imported into sphinx. I think this is a very fun way to write tutorials,
you can almost do a mini study and share it with the world. Every time I
add to the tutorials, I realize there are some new features I would like to
develop.

One thing I would like to try to demonstrate in these tutorials is the ability
of PyBigDFT to work with other code bases. To me, being able to quickly mix
a variety of libraries is the exceptional part of python. In my recent
tutorial on 
[Geometry Optimization](https://l_sim.gitlab.io/bigdft-suite/lessons/GeometryOptimization.html),
I do just this. In the first part, we go
from a SMILES string to an optimized geometry, by integrating 
[openbabel](http://openbabel.org/wiki/Main_Page) into
a BigDFT workflow. In the second part, a slab is created with the aid of
the [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/). In 
[another tutorial](https://l_sim.gitlab.io/bigdft-suite/lessons/MachineLearning.html)
, I integrate PyBigDFT with a machine learning framework called
[qmlcode](https://www.qmlcode.org/index.html). Developing clear datastructures
that enable transformations like this is definitely one of our goals for
PyBigDFT.

> [1] Ratcliff, Laura E., William Dawson, Giuseppe Fisicaro, Damien Caliste,
Stephan Mohr, Augustin Degomme, Brice Videau et al. "Flexibilities of wavelets
as a computational basis set for large-scale electronic structure
calculations." The Journal of Chemical Physics 152, no. 19 (2020): 194110.
