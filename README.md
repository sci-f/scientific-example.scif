# Singularity Example Scientific Workflow

This repository provides a SCI-F (Standard Container Integration Format) Apps example to run a scientific workflow. This analysis was originally used to compare Singularity vs. Docker on cloud providers via [this container](https://github.com/vsoch/singularity-scientific-example) with [interesting results](https://vsoch.github.io/singularity-scientific-example/results/)! If you are interested in a Docker vs. Singularity implementation, see that project. For this small example, we want to give rationale for taking a SCI-F apps approach over a traditional Singularity image. We compare the following equivalent (but different!) implementations:

- [Singularity (without SCI-F)](Singularity.noscif)
- [Singularity (with SCI-F)](Singularity)


## Evaluation
We want to evaluate SCI-F on its ability to expose container metadata, and information about the pipeline and executables inside. First, the goals of SCI-F:

- Containers are **consistent** to allow for comparison. I am able to easily discover relevant software and data for one or more applications defined by the container creator.
- Containers are transparent. If i discover a container and do not have any prior knowledge or metadata, the important executables and metadata are revealed to me.
- Container contents are easily available for **introspection** because they are programmatically parseable. I can run a function over a container, and know exactly the functions available to me, ask for help, and know where to interact with inputs and outputs.
 - Container internal infrastructure is **modular**. Given a set of SCI-F apps from different sources, I can import different install routines and have assurance that environment variables defined for each are sourced correctly for each, and that associated content does not overwrite previous content. Each software and data module must carry, minimally, a unique name and install location in the system.

For each of the above, let's define some basic tests.

### 1. Development Evaluation
For this use case, we are a container developer, and we are writing a singularity build recipe.

**Can I easily define multiple entrypoints**
Singularity standard defaults to a single runscript per container. If I need to expose multiple functions, I either need to write a more complicated single entrypoint to direct the call, or I need to write detailed docs on different exec commands (and executables inside) for the user to run. SCI-F has the advantage here because I don't need to write any logic, the logic is handled by way of creating an "app" for each function.

**Can I easily install files known to modules**
Given that I have two functions for my container to perform, foo and bar, can I generate an install routine that will allow for both shared dependencies (i.e. global) and distinct dependencies? 

A: Singularity standard has one mode to install "global" dependencies, everything goes into the `%post` script and any organization of required data, files, and libraries is scattered around the image. Other than coming up with a manual organization and documenting it, there is no way to cleanly define boundaries that will be easily discovered by the user. With SCI-F, by simply defining an environment, labels, install, or runscript to be in the context of an app, the modularity is automatically generated. When I add a list of files to an app `foo`, I know they are added to the container's predictable location for `foo`. If I add a file to `bin` I know it goes into foo's bin, and is added to the path when `foo` is run. If I add a library to `lib`, I know it is added to `LD_LIBRARY_PATH` when foo is run.

**Can I associate environment and metadata with modules?**
Given two different functions for my container to perform, foo and bar, can I define environment variables and labels (metadata) that I know will be sourced (environment) or exposed (inspect labels) in the context of the app?

A: Singularity standard also has one global shared `%environment`, and `%labels` section. If two functions in the container share the same environment variable and the value is different, this must be resolved manually. With SCI-F, I simply write the different variables to their sections, and have confidence that they will be sourced (environment) or inspected (labels) with clear association to the app.

**Do I need to know standard locations in advance?**
Given that a container has conformance to SCI-F, do I need to know how it works to use it?

A: Instead of requiring the user to know that an app's base might be at `/scif/apps/foo`, we instead expose environment variables (e.g., `SINGULARITY_APPBASE`) that can be referenced at runtime to refer to different apps. This is especially important if, for example, I need to reference the output of one app as input for another, or simply know it's install location.


### 1. Production Evaluation
For this use case, we operate under the scenario that we are familiar with Singularity and the commands to use SCi-F, but we know nothing about the containers. We are a user.

**Do I know what the container does?**
I should be able to run a command to get help or a summary of what the container does (introspection).

A: Yes, for both SCiF and standard Singularity, the creator can define a `%help` section. Both also allow for global `%labels` that might point to documentation, code repositories, or other informative 

**Do I know what executables are important in the container?**
Without much effort, I should have a high level understanding of the different functions that the container performs, as defined by the creator. For example, a container intended for development of variant calling will expose low level tools (e.g, bwa, seqtk) while a container that is intended will expose a pipeline (e.g., mapping).

A: For Singularity standard, if the container performs one function (and one only), then a single runscript / entrypoint is sufficient. Having multiple functions is completely reliant on the programming ability of the creator. If the creator is able to write a more complex runscript that explains different use cases, the container can satisfy this goal. If not, the container is a black box. For SCI-F, without any significant programming ability, the different apps are easily exposed to the user.

**Can I easily get help for an executable?**
Can I ask the container for help, for global help, and for help pertaining to individual functions?

A: For Singularity standard, you **can** ask a container for help, but it's a single global help. For SCI-F apps, you can define help sections for all of the functions that your container serves.


## 1. Build the Image

With Singularity 2.4:

```
sudo singularity build scif.img Singularity
sudo singularity build no-scif.img Singularity.noscif
```
