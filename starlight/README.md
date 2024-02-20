This project packages the starlight project in a container,
this is to aid with running this project in modern
systems. This works with fortran77 which is more than 12 years
old, still the software is useful and commonly used in the astronomy
community.

## Sources

For any reference and docs go to http://www.starlight.ufsc.br/

## How-to

## Build from the repository
```
# Get the source code from GitHub.
git clone https://github.com/ccatalan/astro-containers.git
cd astro-containers

# Go in the starlight folder.
cd starlight

# Build the container.
podman build -t astro-containers/starlight .

#
# If it is required run with -dt instead of -it to
# be executed in detached mode.
#

# Run starlight from this same folder.
podman run --rm -it \
    -v ./spectrum/:/home/starlight/spectrum/:z \
    -v ./mask/:/home/starlight/mask/:z \
    -v ./out/:/home/starlight/out/:z \
    -v ./StCv04.C11.arp220.config:/home/starlight/STARLIGHTv04/StCv04.C11.arp220.config:z \
    localhost/astro-containers/starlight \
        < grid_example.in

# Adjust your data files and execute it as you need.
```

## Reuse the container image

```
# Get the image from quay or github cr
podman pull ghcr.io/ccatalan/astro-containers/starlight:main

# Create the local folders we will mount in the container image
mkdir -p ./spectrum/
mkdir -p ./mask/
mkdir -p ./out/
touch ./StCv04.C11.arp220.config
touch grid_example.in

# For a reference on the content of the files that
# are required to mount see the repository content

#
# If it is required run with -dt instead of -it to 
# be executed in detached mode.
#

# Run the container
podman run --rm -it \
    -v ./spectrum/:/home/starlight/spectrum/:z \
    -v ./mask/:/home/starlight/mask/:z \
    -v ./out/:/home/starlight/out/:z \
    -v ./StCv04.C11.arp220.config:/home/starlight/STARLIGHTv04/StCv04.C11.arp220.config:z \
    ghcr.io/ccatalan/astro-containers/starlight:main \
        < grid_example.in
