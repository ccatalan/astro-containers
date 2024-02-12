This project packages the steckmap project in a container,
this is to aid with running this project in modern
systems.

## Sources

For any reference and docs go to http://astro.u-strasbg.fr/~ocvirk/

Download the steckmap code from http://astro.u-strasbg.fr/~ocvirk/STECKMAP_archive/

The version this container has is STECKMAP_31012011.tar

## How-to

## Build from the repository
```
# Get the source code from GitHub.
git clone https://github.com/ccatalan/astro-containers.git
cd astro-containers

# Go in the starlight folder.
cd steckmap

# Build the container.
podman build -t astro-containers/steckmap .

# Run steckmap example from this same folder.
podman run --rm -it \
    localhost/astro-containers/steckmap \
        -batch STECKMAP/Pierre/POP/sfit.i
```

The previous command should be executed without errors.

## Data

STECKMAP is a tool package aiming at constraining the stellar content and
kinematics from galaxy spectra. To do so, you need two things:

- DATA
- Models against which to compare your data

The following steps are an example about how to map steckmap data and
models to the container.

```
# We create an input folder called 'input'

mkdir input
cd input

# We copy the .i file and the .fits folder to this path
ls
input-> fix_kinematic_BC03_Calzetti_IRAS00188_redshift.i  IRAS00188_redshift
```

Now, the .i file must match the file structure we will mount inside the container

For example the content of the .i file should looks like:

The paths included in the .i are relative to the internal path
structure in the container, not to your local environment, so,
`/home/steckmap/Yorick/STECKMAP` or `STECKMAP` should always be
used. Also the input files inside the input folder are mounted in
`/home/steckmap/input/` so this path always must be also used.

```
include,"sfit.i"
wavel=[4190.0201962423125, 8816.270196242313]
fv="/home/steckmap/input/IRAS00188_redshift/spectrum_0000.fits"
// noplot=1 means no plot!
a=convert_all(fv,noplot=1,z0=0.12834849333333334,SNR0=20)
b=bRbasis3([1.0e5,1.7e10],basisfile="BC03",nbins=30,wavel=wavel,FWHMpix=4.367499800171125)
// np=1 means no plot!
x=sfit(a,b,np=1, kin=0,epar=1,noskip=1,sav=1, RMASK = [[4842.68,4882.68],...,])
print, "##############################################"
print, "Analyzed spectra: 0/191"
print, "##############################################"

quit()
```

Make sure the .i has for convert_all and sfit, noplot=1, and np=1
respectively as in the example.
All .i files MUST end with `quit()` as the container
entrypoint will be calling Yorick so we need the
container to end once the input files are executed.

From the parent folder of the input path (the steckmap folder) execute:



## Running as a detached service (no additional x11 forward is required)

We have all the requirements locally, so no need to extra configuration parameters.

```
podman run --rm -it \
    -v ./input/:/home/steckmap/input/:z \
    localhost/astro-containers/steckmap \
        -batch /home/steckmap/input/fix_kinematic_BC03_Calzetti_IRAS00188_redshift.i
```

This will mount the input folder in the place where the .i file expects to find all the
other sources all inside the container.

```
# Get the image from quay or github cr
podman pull quay.io/ccatalan/astro-containers-steckmap:latest
```

## OLD - DEPRECATED without the X11 server:

It is needed to `-v /tmp/.X11-unix:/tmp/.X11-unix` to mount the X11 socket,
also `-e DISPLAY=unix$DISPLAY` to  pass the display,
and `--security-opt label=type:container_runtime_t` to allow run.

```
podman run --rm -it \
    -e DISPLAY=unix$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --security-opt label=type:container_runtime_t  \
    -v ./input/:/home/steckmap/input/:z \
    localhost/astro-containers/steckmap \
        -batch /home/steckmap/input/fix_kinematic_BC03_Calzetti_IRAS00188_redshift.i
```

