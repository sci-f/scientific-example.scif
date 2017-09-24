# Singularity Example Scientific Workflow

This repository provides a SCI-F (Standard Container Integration Format) Apps example to run a scientific workflow. This analysis was originally used to compare Singularity vs. Docker on cloud providers via [this container](https://github.com/vsoch/singularity-scientific-example) with [interesting results](https://vsoch.github.io/singularity-scientific-example/results/)! If you are interested in a Docker vs. Singularity implementation, see that project. For this small example, we want to give rationale for taking a SCI-F apps approach over a traditional Singularity image. We compare the following equivalent (but different!) implementations:

- [Singularity (without SCI-F)](Singularity.noscif)
- [Singularity (with SCI-F)](Singularity)
