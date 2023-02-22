# pip version >= 10
from setuptools import setup, find_packages


with open("requirements.txt", "r") as f:
    install_requires = f.read().split("\n")

# Read __version__ variable
__version__ = "0.0.0"
exec(open("topology_loop_caller/__init__.py").read())

setup(
    name="topology_loop_caller",
    version=__version__,
    packages=find_packages(),
    url="",
    license="",
    author="Alexander_Kuznetsov",
    author_email="alexander.kuznetsov.bioinf@gmail.com",
    description="Persistent homology - based pipeline of feature generation and loop calling",
    install_requires=install_requires,
    package_dir={"topology_loop_caller": "topology_loop_caller"},
    include_package_data=True,
)
