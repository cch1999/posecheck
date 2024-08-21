from setuptools import setup, find_packages

# Read requirements.txt and store its contents in a list
with open('requirements.txt') as f:
    required = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="posecheck",
    version="1.3.1",
    description="A library for benchmarking poses of 3D SBDD models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Charles Harris",
    author_email="cch57@cam.ac.uk",
    packages=find_packages(include=["posecheck",
                                    "posecheck.*"]),
    install_requires=required,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
