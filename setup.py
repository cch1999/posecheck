from setuptools import setup, find_packages

setup(
    name="posecheck",
    version="1.1",
    description="A library for pose estimation benchmarks",
    author="Charles Harris",
    author_email="cch57@cam.ac.uk",
    packages=find_packages(include=["posecheck", "posecheck.*"]),
    # install_requires=[
    #    "numpy>=1.19.5",
    #    "tensorflow>=2.6",
    #    # Add other dependencies here
    # ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
