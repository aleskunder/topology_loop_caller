# pip version >= 10
from setuptools import setup, find_packages


with open("requirements.txt", "r") as f:
    install_requires = f.read().split("\n")

with open("README.md", "r") as fh:
    long_description = fh.read()

# Read __version__ variable
__version__ = "0.0.0"
exec(open("topology_loop_caller/__init__.py").read())

setup(
    name="topology_loop_caller",
    version=__version__,
    packages=find_packages(),
    url="https://github.com/aleskunder/topology_loop_caller",
    license="",  # Add after clarif
    author="Khrameeva Lab",  # Clarif
    author_email="alexander.kuznetsov.bioinf@gmail.com",
    description="Persistent homology - based pipeline of feature generation and loop calling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=install_requires,
    package_dir={"topology_loop_caller": "topology_loop_caller"},
    include_package_data=True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
