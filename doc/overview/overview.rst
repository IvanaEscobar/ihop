Overview
========

This document seves as a guide for basic usage and theory of underwater 
acoustics within a general ocean circulation model, MITgcm. It gives a 
description of a simple use case via an example of a forward in time run.

Introduction
============

**ihop** is a simulation of underwater acoustics based on ray-traced wavefronts.
The code started from BELLHOP, by Acoustics Toolbox :cite:`MP11`, and has been retrofit to
work within MITgcm :cite:`AA18`. 
